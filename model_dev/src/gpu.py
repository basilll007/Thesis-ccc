"""GPU implementations of the Phase 2 baselines (scaler, ridge, gradient boosting) with small host-RAM footprint.

What: feature blocks live on the GPU as float32 tensors; standardisation, closed-form ridge (all alphas from one
eigendecomposition) and the MLP run on the GPU; gradient boosting uses XGBoost's CUDA backend.

Why (Phase 2 plan, Amendment 1): the CPU runs held several GB of feature copies in host RAM and were stopped by the
host for low system memory. Moving the matrices to the 8 GB GPU keeps host RAM small and follows the instruction to
train on the GPU. Ridge here is numerically the same estimator as sklearn `Ridge(alpha, fit_intercept=True)`
(centred least squares with an L2 penalty), solved in float64; `tests/test_gpu.py` checks the agreement.
"""
from __future__ import annotations

import numpy as np
import torch


def device() -> str:
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required for model_dev training (GPU-only policy); no CPU fallback for training")
    return "cuda"


def to_gpu(a: np.ndarray, dtype=torch.float32) -> torch.Tensor:
    return torch.as_tensor(np.ascontiguousarray(a), dtype=dtype, device=device())


def fit_scaler(X: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    """Column mean/sd over the given (training) rows, ignoring NaN; constant columns get sd 1."""
    mask = ~torch.isnan(X)
    n = mask.sum(0).clamp(min=1)
    Xz = torch.where(mask, X, torch.zeros((), device=X.device, dtype=X.dtype))
    mu = Xz.sum(0) / n
    var = (torch.where(mask, (X - mu) ** 2, torch.zeros((), device=X.device, dtype=X.dtype))).sum(0) / n
    sd = var.sqrt()
    sd = torch.where((sd == 0) | ~torch.isfinite(sd), torch.ones_like(sd), sd)
    return mu, sd


def apply_scaler(X: torch.Tensor, mu: torch.Tensor, sd: torch.Tensor) -> torch.Tensor:
    """Standardise; NaN (empty tile) becomes 0, i.e. the training mean."""
    return torch.nan_to_num((X - mu) / sd, nan=0.0)


class RidgePath:
    """Ridge for many alphas at once: centred normal equations, one eigendecomposition, float64 on the GPU."""

    def __init__(self, X: torch.Tensor, Y: torch.Tensor, chunk: int = 8192):
        # float64 Gram matrix accumulated over row chunks: no full-size float64 copy of X is ever materialised
        self.xm = X.double().mean(0) if len(X) <= chunk else sum(X[i:i + chunk].double().sum(0) for i in range(0, len(X), chunk)) / len(X)
        self.ym = Y.double().mean(0)
        p, k = X.shape[1], Y.shape[1]
        G = torch.zeros((p, p), dtype=torch.float64, device=X.device)
        XtY = torch.zeros((p, k), dtype=torch.float64, device=X.device)
        for i in range(0, len(X), chunk):
            xc = X[i:i + chunk].double() - self.xm
            G += xc.T @ xc
            XtY += xc.T @ (Y[i:i + chunk].double() - self.ym)
            del xc
        self.evals, self.V = torch.linalg.eigh(G)
        del G
        self.b = self.V.T @ XtY                              # [p, k]

    def coef(self, alpha: float) -> torch.Tensor:
        return self.V @ (self.b / (self.evals.clamp(min=0) + alpha)[:, None])

    def predict(self, X: torch.Tensor, alpha: float, chunk: int = 16384) -> torch.Tensor:
        W = self.coef(alpha)
        return torch.cat([((X[i:i + chunk].double() - self.xm) @ W + self.ym).float() for i in range(0, len(X), chunk)])


def xgb_fit_select(Xtr: np.ndarray, ytr: np.ndarray, Xval: np.ndarray, yval: np.ndarray, lr: float, leaves: int,
                   max_rounds: int, seed: int):
    """XGBoost (CUDA) with the round count chosen on the inner-validation RMSE; returns (best_rounds, val_mse)."""
    import xgboost as xgb
    params = {"tree_method": "hist", "device": "cuda", "learning_rate": lr, "max_leaves": leaves, "max_depth": 0,
              "grow_policy": "lossguide", "min_child_weight": 50, "objective": "reg:squarederror", "seed": seed,
              "eval_metric": "rmse"}
    dtr, dval = xgb.DMatrix(Xtr, label=ytr), xgb.DMatrix(Xval, label=yval)
    ev: dict = {}
    xgb.train(params, dtr, max_rounds, evals=[(dval, "val")], evals_result=ev, verbose_eval=False)
    rmse = np.array(ev["val"]["rmse"])
    best = int(np.argmin(rmse)) + 1
    return best, float(rmse[best - 1] ** 2)


def xgb_fit_predict(Xtr: np.ndarray, ytr: np.ndarray, Xte: np.ndarray, lr: float, leaves: int, rounds: int, seed: int) -> np.ndarray:
    import xgboost as xgb
    params = {"tree_method": "hist", "device": "cuda", "learning_rate": lr, "max_leaves": leaves, "max_depth": 0,
              "grow_policy": "lossguide", "min_child_weight": 50, "objective": "reg:squarederror", "seed": seed}
    booster = xgb.train(params, xgb.DMatrix(Xtr, label=ytr), rounds)
    return booster.predict(xgb.DMatrix(Xte))


class GPUMLP:
    """Multi-output MLP with early stopping on the inner-validation fold; inputs already on the GPU."""

    def __init__(self, n_in: int, n_out: int, hidden=(256, 128), dropout=0.1, lr=1e-3, wd=1e-4, seed=0, batch=1024,
                 max_epochs=100, patience=10):
        import torch.nn as nn
        torch.manual_seed(seed); torch.cuda.manual_seed_all(seed)
        layers, d = [], n_in
        for h in hidden:
            layers += [nn.Linear(d, h), nn.ReLU(), nn.Dropout(dropout)]
            d = h
        layers.append(nn.Linear(d, n_out))
        self.net = nn.Sequential(*layers).to(device())
        self.lr, self.wd, self.batch, self.max_epochs, self.patience, self.seed = lr, wd, batch, max_epochs, patience, seed
        self.history: list[dict] = []

    def fit(self, Xtr: torch.Tensor, ytr: torch.Tensor, Xval: torch.Tensor, yval: torch.Tensor):
        opt = torch.optim.AdamW(self.net.parameters(), lr=self.lr, weight_decay=self.wd)
        g = torch.Generator(device="cuda").manual_seed(self.seed)
        best, bad, best_state = np.inf, 0, None
        for epoch in range(self.max_epochs):
            self.net.train()
            perm = torch.randperm(len(Xtr), generator=g, device="cuda")
            tot = torch.zeros((), device="cuda")
            for i in range(0, len(perm), self.batch):
                idx = perm[i:i + self.batch]
                opt.zero_grad(set_to_none=True)
                loss = ((self.net(Xtr[idx]) - ytr[idx]) ** 2).mean()
                loss.backward()
                opt.step()
                tot += loss.detach() * len(idx)
            vl = float(((self.predict(Xval) - yval) ** 2).mean())
            self.history.append({"epoch": epoch, "train_mse": float(tot) / len(Xtr), "val_mse": vl})
            if vl < best - 1e-6:
                best, bad = vl, 0
                best_state = {k: v.detach().clone() for k, v in self.net.state_dict().items()}
            else:
                bad += 1
                if bad >= self.patience:
                    break
        self.net.load_state_dict(best_state)
        self.best_val, self.epochs = best, len(self.history)
        return self

    def fit_fixed(self, Xtr: torch.Tensor, ytr: torch.Tensor, epochs: int):
        """Fixed epoch budget (used for the external refit, where no validation fold exists)."""
        opt = torch.optim.AdamW(self.net.parameters(), lr=self.lr, weight_decay=self.wd)
        g = torch.Generator(device="cuda").manual_seed(self.seed)
        for _ in range(epochs):
            self.net.train()
            perm = torch.randperm(len(Xtr), generator=g, device="cuda")
            for i in range(0, len(perm), self.batch):
                idx = perm[i:i + self.batch]
                opt.zero_grad(set_to_none=True)
                ((self.net(Xtr[idx]) - ytr[idx]) ** 2).mean().backward()
                opt.step()
        self.epochs = epochs
        return self

    @torch.no_grad()
    def predict(self, X: torch.Tensor) -> torch.Tensor:
        self.net.eval()
        return torch.cat([self.net(X[i:i + 16384]) for i in range(0, len(X), 16384)])
