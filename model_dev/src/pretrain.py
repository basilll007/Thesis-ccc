"""Phase 6: masked-expression pretraining feasibility (exploratory, development set only, fold-disciplined).

What (plan: model_dev/plans/phase6_pretraining.md): inside each outer fold, pretrain a padded-tree encoder on the
inner-training cells of all types to predict 20 % masked genes of a cell from its remaining genes, covariates and
2-hop neighbourhood; early stop on the inner-validation masked loss. Then fit a ridge probe on the frozen seed
embedding (target genes masked; a second 'neighbourhood-only' probe masks all own genes) for the Phase 1 targets,
and compare with the Phase 2/3 models. Also a learning curve over the number of training tumours.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402
import assemble as A  # noqa: E402
import features as F  # noqa: E402
import models as M  # noqa: E402
import train as T  # noqa: E402
from baselines import fit_ridge, merge_pairs, pooled_r2, sse_sst_by_group  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
import torch.nn as nn  # noqa: E402

PHASE = "phase6"


class MaskedEncoder(nn.Module):
    """Tree encoder whose seed input is the cell's own (partially masked) expression plus covariates."""

    def __init__(self, n_genes: int, n_classes: int, hidden: int, emb: int, dropout: float):
        super().__init__()
        self.tree = M.TreeModel(n_genes, n_classes, n_genes + 2 + n_genes, hidden, hidden, emb, "mean", dropout)
        self.tree.head = nn.Identity()                                   # expose the hidden state
        self.decoder = nn.Linear(hidden, n_genes)

    def embed(self, b: M.TreeBatch) -> torch.Tensor:
        return self.tree(b)

    def forward(self, b: M.TreeBatch) -> torch.Tensor:
        return self.decoder(self.embed(b))


def masked_inputs(x: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    """Own-expression input: masked genes set to 0, plus the mask indicator (so the model knows what is hidden)."""
    return torch.cat([x * (~mask), mask.float()], -1)


def pretrain_fold(store, seeds_tr: np.ndarray, seeds_val: np.ndarray, cov_all: torch.Tensor, cfg: dict, n_classes: int,
                  device: str, seed: int, log):
    C.set_seed(seed)
    enc = MaskedEncoder(store.X.shape[1], n_classes, cfg["hidden"], cfg["emb"], cfg["dropout"]).to(device)
    opt = torch.optim.AdamW(enc.parameters(), lr=cfg["lr"], weight_decay=cfg["weight_decay"])
    gen = torch.Generator(device=device).manual_seed(seed)
    k1, k2, B = cfg["k1"], cfg["k2"], cfg["batch"]
    hist, best, bad, best_state = [], np.inf, 0, None

    def step(idx: np.ndarray, train: bool) -> float:
        enc.train(train)
        tot, n = 0.0, 0
        with torch.set_grad_enabled(train):
            for i in range(0, len(idx), B):
                g = torch.as_tensor(idx[i:i + B]).long().to(device)
                x = store.X[g]
                mask = torch.rand(x.shape, generator=gen, device=device) < cfg["mask_frac"]
                seed_cov = torch.cat([cov_all[g], masked_inputs(x, mask)], -1)
                tb = M.TreeBatch(store, g, seed_cov, -1, "full", k1, k2)
                with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=device == "cuda"):
                    pred = enc(tb).float()
                loss = (((pred - x) ** 2) * mask).sum() / mask.sum().clamp(min=1)
                if train:
                    opt.zero_grad(set_to_none=True); loss.backward(); opt.step()
                tot += float(loss) * len(g); n += len(g)
        return tot / max(n, 1)

    rng = np.random.default_rng(seed)
    for epoch in range(cfg["max_epochs"]):
        tr = rng.permutation(seeds_tr)[:cfg["cells_per_epoch"]]
        tl = step(tr, True)
        vl = step(seeds_val[:cfg["val_cells"]], False)
        hist.append({"epoch": epoch, "train_masked_mse": tl, "val_masked_mse": vl})
        log.info("pretrain epoch %d: train %.4f val %.4f", epoch, tl, vl)
        if vl < best - 1e-6:
            best, bad, best_state = vl, 0, {k: v.detach().clone() for k, v in enc.state_dict().items()}
        else:
            bad += 1
            if bad >= cfg["patience"]:
                break
    enc.load_state_dict(best_state)
    return enc, hist


def embed_cells(enc, store, g: torch.Tensor, cov: torch.Tensor, target_gene_idx: np.ndarray, mask_all: bool, cfg, device):
    enc.eval()
    out = []
    gmask = torch.zeros(store.X.shape[1], dtype=torch.bool, device=device)
    gmask[torch.as_tensor(target_gene_idx).to(device)] = True
    with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=device == "cuda"):
        for i in range(0, len(g), 2048):
            gi = g[i:i + 2048]
            x = store.X[gi]
            mask = torch.ones_like(x, dtype=torch.bool) if mask_all else gmask[None, :].expand_as(x)
            tb = M.TreeBatch(store, gi, torch.cat([cov[gi], masked_inputs(x, mask)], -1), -1, "full", cfg["k1"], cfg["k2"])
            out.append(enc.embed(tb).float().cpu().numpy())
    return np.vstack(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(C.MODEL_DEV / "configs" / "pretrain.yaml"))
    args = ap.parse_args()
    cfg = C.load_config(args.config)
    log = C.start_run(PHASE, "pretrain", cfg)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    fold_table = A.folds("grouped"); n_folds = int(fold_table.fold.max()) + 1
    spec = A.targets_spec()
    store0, offsets, classes = T.build_store("development", max(cfg["k1"], cfg["k2"]), device, log)
    genes = [str(g) for g in A.load_section(A.section_files("development")[0])["genes"]]
    rows, pairs = [], {}
    curves = []
    for k in range(n_folds):
        t0 = time.time()
        secs_tr = set(fold_table[(fold_table.fold != k) & (fold_table.fold != (k + 1) % n_folds)].section)
        secs_val = set(fold_table[fold_table.fold == (k + 1) % n_folds].section)
        secs_all_tr = secs_tr | secs_val
        gmask = torch.as_tensor(np.isin(store0.section_of, list(secs_all_tr))).to(device)
        mu, sd, cmu, csd = T.fold_standardisation(store0, gmask)
        store = T.StandardisedStore(store0, mu, sd, cmu, csd)
        cov_all = store.cov[torch.arange(store0.n, device=device)]
        idx_tr = np.where(np.isin(store0.section_of, list(secs_tr)))[0]
        idx_val = np.where(np.isin(store0.section_of, list(secs_val)))[0]
        enc, hist = pretrain_fold(store, idx_tr, np.random.default_rng(k).permutation(idx_val), cov_all, cfg, len(classes), device, cfg["seed"], log)
        curves += [{**h, "fold": k} for h in hist]
        ck = C.CACHE / "checkpoints" / PHASE / f"encoder_fold{k}.pt"; ck.parent.mkdir(parents=True, exist_ok=True)
        torch.save({"state": enc.state_dict(), "mu": mu.cpu(), "sd": sd.cpu(), "cfg": cfg}, ck)
        for receiver in cfg["receivers"]:
            s = spec[receiver]
            data = A.build_receiver("development", receiver, log)
            tr_in, val, tr_all, test = A.fold_masks(data, fold_table, k, n_folds)
            g = torch.as_tensor(np.array([offsets[x] + r for x, r in zip(data.section, data.row)])).long().to(device)
            tidx = np.array([genes.index(x) for x in s["genes"]])
            ytr_mean = data.Y[tr_all].mean(axis=0)
            for probe in ("target_masked", "neighbourhood_only"):
                E = embed_cells(enc, store, g, cov_all, tidx, probe == "neighbourhood_only", cfg, device)
                fs = F.FeatureSet(rows=np.arange(len(E)))
                fs.add("emb", E, [f"emb_{i}" for i in range(E.shape[1])], [F.Provenance("tile")] * E.shape[1])  # other cells' + own non-target genes
                emu, esd = F.fit_scaler(E[tr_in]); Ein, Ev = F.apply_scaler(E[tr_in], emu, esd), F.apply_scaler(E[val], emu, esd)
                vs = {a: float(((fit_ridge(Ein, data.Y[tr_in], a).predict(Ev) - data.Y[val]) ** 2).mean()) for a in cfg["ridge_alphas"]}
                a = min(vs, key=vs.get)
                emu2, esd2 = F.fit_scaler(E[tr_all])
                pred = fit_ridge(F.apply_scaler(E[tr_all], emu2, esd2), data.Y[tr_all], a).predict(F.apply_scaler(E[test], emu2, esd2))
                n_prog = len(s["programs"])
                for j, t in enumerate(data.target_names + ["genes_pooled"]):
                    if t == "genes_pooled":
                        yj, pj, mj = data.Y[test][:, n_prog:], pred[:, n_prog:], ytr_mean[n_prog:]
                    else:
                        yj, pj, mj = data.Y[test][:, j:j + 1], pred[:, j:j + 1], ytr_mean[j:j + 1]
                    p = sse_sst_by_group(yj, pj, data.group[test], mj)
                    key = (receiver, probe, t)
                    pairs[key] = merge_pairs([pairs.get(key, {}), p])
                    rows.append({"receiver": receiver, "probe": probe, "fold": k, "target": t, "alpha": a, "test_r2_fold": pooled_r2(p)})
        log.info("fold %d done in %.0f s", k, time.time() - t0)
    res = pd.DataFrame(rows)
    C.save_table(res, PHASE, "probe_results.csv", log)
    pd.DataFrame(curves).to_csv(C.RESULTS / PHASE / "pretrain_curves.csv", index=False)
    # comparisons
    summ2 = pd.read_csv(C.RESULTS / "phase2" / "baseline_summary.csv")
    g3 = pd.read_csv(C.RESULTS / "phase3" / "graph_summary.csv")
    out = []
    for (receiver, probe, t), p in pairs.items():
        r2, lo, hi = C.tumour_bootstrap(p, cfg["n_boot"])
        bname, bp = T.baseline_pairs(receiver, t)
        gb, gblo, gbhi = C.tumour_bootstrap_diff(bp, p, cfg["n_boot"])
        gp = T.pairs_from_csv("phase3", "graph", [receiver])
        sel = g3[(g3.receiver == receiver) & (g3.target == t)].sort_values("mean_inner_val_mse").iloc[0]
        gpairs = T.seed_mean_pairs(gp, receiver, "full", sel.arch, float(sel.selected_lr), [0, 1, 2], t)
        gg, gglo, gghi = C.tumour_bootstrap_diff(gpairs, p, cfg["n_boot"])
        out.append({"receiver": receiver, "probe": probe, "target": t, "r2": r2, "ci_low": lo, "ci_high": hi,
                    "strongest_baseline": bname, "probe_minus_baseline": gb, "pb_ci_low": gblo, "pb_ci_high": gbhi,
                    "graph_model": sel.arch, "probe_minus_graph": gg, "pg_ci_low": gglo, "pg_ci_high": gghi,
                    "pretraining_helps": bool(probe == "neighbourhood_only" and gglo > 0)})
    summ = pd.DataFrame(out)
    C.save_table(summ, PHASE, "probe_summary.csv", log)
    log.info("probe summary:\n%s", summ.round(4).to_string(index=False))
    log.info("done")


if __name__ == "__main__":
    main()
