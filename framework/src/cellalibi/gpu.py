"""Optional GPU (PyTorch/CUDA) backend for the permutation nulls in `cellalibi.contact`.

Used when `AuditConfig(device="cuda")`. The statistics, null designs and add-one p-values are identical to the CPU
path; only the random-number stream differs (torch generator instead of numpy), so p-values agree up to Monte-Carlo
error. Grouped (within-tile / within-type) permutations use float64 sort keys = group code + U(0, 1), exactly as the
CPU path does, so every within-group ordering is equally likely.
"""
from __future__ import annotations

import numpy as np

try:  # torch is optional
    import torch
except Exception:  # pragma: no cover
    torch = None


def available() -> bool:
    return torch is not None and torch.cuda.is_available()


class GpuPerm:
    """Holds a seeded torch generator on the device and provides batched permutation nulls."""

    def __init__(self, seed: int, device: str = "cuda", chunk: int = 250):
        if not available():
            raise RuntimeError("CUDA not available")
        self.dev = torch.device(device)
        self.gen = torch.Generator(device=self.dev)
        self.gen.manual_seed(int(seed) & 0x7FFFFFFFFFFFFFFF)
        self.chunk = chunk

    def t(self, x, dtype=None):
        return torch.as_tensor(np.asarray(x), device=self.dev, dtype=dtype)

    def _grouped_perm(self, groups_t, b):
        n = groups_t.shape[0]
        key = groups_t.to(torch.float64)[None, :] + torch.rand((b, n), generator=self.gen, device=self.dev, dtype=torch.float64)
        order = torch.argsort(key, dim=1)
        pos = torch.argsort(groups_t, stable=True)
        out = torch.empty_like(order)
        out[:, pos] = order
        return out

    def contact_nulls(self, Ls, Rr, es, er, gs, gr, n_perms):
        Ls, Rr = self.t(Ls, torch.float64), self.t(Rr, torch.float64)
        es, er = self.t(es, torch.int64), self.t(er, torch.int64)
        gs, gr = self.t(gs, torch.int64), self.t(gr, torch.int64)
        out = []
        done = 0
        while done < n_perms:
            b = min(self.chunk, n_perms - done)
            ps = self._grouped_perm(gs, b)
            pr = self._grouped_perm(gr, b)
            Lp = torch.gather(Ls.expand(b, -1), 1, ps)
            Rp = torch.gather(Rr.expand(b, -1), 1, pr)
            out.append((Lp[:, es] * Rp[:, er]).mean(1))
            done += b
        return torch.cat(out).cpu().numpy()

    def ring_null(self, vals, is_x, groups, n_perms):
        v = self.t(vals, torch.float64)
        x = self.t(is_x, torch.bool)
        g = self.t(groups, torch.int64)
        nx, nc = x.sum(), (~x).sum()
        out = []
        done = 0
        while done < n_perms:
            b = min(self.chunk, n_perms - done)
            p = self._grouped_perm(g, b)
            V = torch.gather(v.expand(b, -1), 1, p)
            out.append(V[:, x].sum(1) / nx - V[:, ~x].sum(1) / nc)
            done += b
        return torch.cat(out).cpu().numpy()

    def label_matched_null(self, Llog, Rlog, n_s, n_r, n_perms):
        L, R = self.t(Llog, torch.float64), self.t(Rlog, torch.float64)
        n = L.shape[0]
        out = []
        done = 0
        b0 = max(1, min(self.chunk, int(2e7 // max(n, 1))))
        while done < n_perms:
            b = min(b0, n_perms - done)
            idx = torch.argsort(torch.rand((b, n), generator=self.gen, device=self.dev), dim=1)
            ls = torch.gather(L.expand(b, -1), 1, idx[:, :n_s]).mean(1)
            rs = torch.gather(R.expand(b, -1), 1, idx[:, n_s:n_s + n_r]).mean(1)
            out.append(ls * rs)
            done += b
        return torch.cat(out).cpu().numpy()

    def phi_null(self, is_t, Mk, den, n_perms):
        it = self.t(is_t, torch.bool)
        M = self.t(Mk, torch.float64)
        d = self.t(den, torch.float64)
        n = it.shape[0]
        out = []
        done = 0
        while done < n_perms:
            b = min(self.chunk, n_perms - done)
            keys = torch.argsort(torch.rand((b, n), generator=self.gen, device=self.dev), dim=1)
            mt = it[keys].to(torch.float64)
            nt = mt.sum(1, keepdim=True)
            ni = (1 - mt).sum(1, keepdim=True)
            m_t = (mt @ M) / nt
            m_i = ((1 - mt) @ M) / ni
            ratio = (m_t - m_i) / torch.clamp(d[None, :] - m_i, min=1e-6)
            out.append(torch.quantile(ratio, 0.5, dim=1))
            done += b
        return torch.cat(out).cpu().numpy()
