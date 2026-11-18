"""Claim-level contact audit (v2 engine).

What this module does:
For one claim "sender type S signals to receiver type R through ligand L and receptor Rg" it computes,
on a symmetric k-NN contact graph:
  * contact test (global): Q = E/M; null permutes L among S cells and Rg among R cells, so cell-type
    abundance and type-level expression are preserved and only who-touches-whom is broken;
  * contact test (block): the same permutation restricted to square tiles (spatially constrained null);
  * ring-control test: Rg in receivers touching an L+ sender versus nearby receivers that do not touch one;
  * naive abundance score: matched baseline M under global cell-type label permutation (thesis style);
  * plant-and-shift power for a declared effect size;
  * claim-specific segmentation-spillover check (sender-marker contamination of touching receivers and
    whether it can account for the observed receptor excess).

Why:
v1 reported only the first failing check, ran the spillover diagnostic for one pair and applied it to all
claims, and used a label-permutation null whose reference distribution mixes expression across cell types.
Plan: framework/plans/v2_phase2.md.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np
import scipy.sparse as sp
from sklearn.neighbors import KDTree, NearestNeighbors


@dataclass
class AuditConfig:
    k: int = 6
    n_perms: int = 1000
    alpha: float = 0.05
    tile_um: float = 100.0
    ring_factor: float = 4.0
    power_effect: float = 1.5
    power_floor: float = 0.5
    min_receivers: int = 25
    min_edges: int = 10
    min_ring_group: int = 5
    spill_marker_ratio: float = 10.0
    spill_marker_min_cp10k: float = 1.0
    spill_n_markers: int = 10
    spill_min_markers: int = 3
    spill_min_isolated: int = 10
    spill_phi_min: float = 0.005
    spill_explain_frac: float = 0.5
    perm_chunk: int = 100
    seed: int = 42

    @property
    def min_p(self) -> float:
        return 1.0 / (self.n_perms + 1)


class SpatialContext:
    """Geometry shared by every claim in a section: contact graph, diameter proxy, tiles."""

    def __init__(self, xy: np.ndarray, labels: Sequence[str], cfg: AuditConfig | None = None):
        self.cfg = cfg or AuditConfig()
        self.xy = np.asarray(xy, dtype=float)
        self.labels = np.asarray(labels).astype(str)
        n = self.xy.shape[0]
        self.n = n
        k = min(self.cfg.k, max(n - 1, 1))
        nn = NearestNeighbors(n_neighbors=k + 1).fit(self.xy)
        dist, idx = nn.kneighbors(self.xy)
        a = np.repeat(np.arange(n), k)
        b = idx[:, 1:].ravel()
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        key = np.unique(lo.astype(np.int64) * n + hi)
        self.pa = (key // n).astype(np.int64)
        self.pb = (key % n).astype(np.int64)
        nd = dist[:, 1]
        self.diameter = float(np.median(nd[nd > 0])) if np.any(nd > 0) else 10.0
        rng = np.random.default_rng(self.cfg.seed)
        off = rng.uniform(0, self.cfg.tile_um, size=2)
        g = np.floor((self.xy - self.xy.min(0) + off) / self.cfg.tile_um).astype(np.int64)
        _, self.tiles = np.unique(g[:, 0] * (g[:, 1].max() + 1) + g[:, 1], return_inverse=True)
        self.tiles = self.tiles.astype(np.int64)

    def directed_pairs(self, s_mask: np.ndarray, r_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Contact pairs (sender cell, receiver cell) from the symmetric graph."""
        pa, pb = self.pa, self.pb
        f1 = s_mask[pa] & r_mask[pb]
        f2 = s_mask[pb] & r_mask[pa]
        return np.concatenate([pa[f1], pb[f2]]), np.concatenate([pb[f1], pa[f2]])


class ExpressionData:
    """Counts -> CP10K (lin) and log1p(CP10K) (log); per-gene dense access; full lin for markers."""

    def __init__(self, counts, gene_names: Sequence[str], is_log1p_cp10k: bool = False):
        self.genes = [str(g) for g in gene_names]
        self.gi = {g: i for i, g in enumerate(self.genes)}
        M = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(np.asarray(counts))
        M = M.astype(np.float64)
        if is_log1p_cp10k:
            lin = M.copy()
            lin.data = np.expm1(lin.data)
            self.counts = None
        else:
            self.counts = M
            tot = np.asarray(M.sum(axis=1)).ravel()
            scale = np.where(tot > 0, 1e4 / np.maximum(tot, 1e-12), 0.0)
            lin = (sp.diags(scale) @ M).tocsr()
        self.lin = lin
        self._cache: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    def gene(self, g: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (lin, log, positive) dense vectors for gene g."""
        if g not in self._cache:
            j = self.gi[g]
            lin = np.asarray(self.lin[:, j].todense()).ravel()
            if self.counts is not None:
                pos = np.asarray(self.counts[:, j].todense()).ravel() > 0
            else:
                pos = lin > 0
            self._cache[g] = (lin, np.log1p(lin), pos)
        return self._cache[g]


def _grouped_perm_keys(groups: np.ndarray, rng: np.random.Generator, b: int) -> np.ndarray:
    """(b, n) index arrays that permute positions within groups; rows map sorted-group slots to elements."""
    n = groups.shape[0]
    key = groups[None, :].astype(np.float64) + rng.random((b, n))
    order = np.argsort(key, axis=1, kind="stable")
    pos = np.argsort(groups, kind="stable")
    out = np.empty((b, n), dtype=np.int64)
    out[:, pos] = order
    return out


def _addone_p(null: np.ndarray, obs: float) -> float:
    if not np.isfinite(obs):
        return float("nan")
    return float((np.sum(null >= obs - 1e-12) + 1) / (null.shape[0] + 1))


def _shift_power(null: np.ndarray, shift: float, alpha: float) -> float:
    crit = np.quantile(null, 1 - alpha)
    return float(np.mean(null + shift > crit))


def contact_nulls(Ls: np.ndarray, Rr: np.ndarray, es: np.ndarray, er: np.ndarray, gs: np.ndarray, gr: np.ndarray,
                  n_perms: int, rng: np.random.Generator, chunk: int) -> np.ndarray:
    """Null E values: L permuted within groups gs over senders, Rg within groups gr over receivers."""
    out = np.empty(n_perms)
    done = 0
    while done < n_perms:
        b = min(chunk, n_perms - done)
        ps = _grouped_perm_keys(gs, rng, b)
        pr = _grouped_perm_keys(gr, rng, b)
        Lp = Ls[ps]
        Rp = Rr[pr]
        out[done:done + b] = np.mean(Lp[:, es] * Rp[:, er], axis=1)
        done += b
    return out


def ring_null(vals: np.ndarray, is_x: np.ndarray, groups: np.ndarray, n_perms: int, rng: np.random.Generator,
              chunk: int) -> np.ndarray:
    out = np.empty(n_perms)
    nx = is_x.sum()
    nc = (~is_x).sum()
    done = 0
    while done < n_perms:
        b = min(chunk, n_perms - done)
        p = _grouped_perm_keys(groups, rng, b)
        V = vals[p]
        out[done:done + b] = V[:, is_x].sum(1) / nx - V[:, ~is_x].sum(1) / nc
        done += b
    return out


def label_matched_null(Llog: np.ndarray, Rlog: np.ndarray, n_s: int, n_r: int, n_perms: int,
                       rng: np.random.Generator) -> np.ndarray:
    """Thesis-style matched baseline M under global cell-type label permutation."""
    n = Llog.shape[0]
    out = np.empty(n_perms)
    for k in range(n_perms):
        p = rng.permutation(n)
        out[k] = Llog[p[:n_s]].mean() * Rlog[p[n_s:n_s + n_r]].mean()
    return out


def spillover_check(ctx: SpatialContext, ex: ExpressionData, s_mask: np.ndarray, r_mask: np.ndarray,
                    ligand: str, receptor: str, cfg: AuditConfig, rng: np.random.Generator) -> dict[str, Any]:
    """Claim-specific spillover diagnostic (plan v2_phase2 section 2.2)."""
    out: dict[str, Any] = {"spill_status": "NOT_RUN"}
    si, rj = ctx.directed_pairs(s_mask, r_mask)
    touching = np.zeros(ctx.n, bool)
    touching[rj] = True
    rec = np.flatnonzero(r_mask)
    T = rec[touching[rec]]
    if s_mask.sum() == 0 or len(T) < 5:
        out["spill_reason"] = "too few touching receivers"
        return out
    tree = KDTree(ctx.xy[s_mask])
    dS = tree.query(ctx.xy[rec], k=1)[0].ravel()
    nonT = ~touching[rec]
    I = rec[nonT & (dS > 2 * ctx.diameter)]
    out["spill_isolated_rule"] = "no S contact and nearest S > 2d"
    if len(I) < cfg.spill_min_isolated:
        I = rec[nonT]
        out["spill_isolated_rule"] = "no S contact (fallback)"
    if len(I) < cfg.spill_min_isolated:
        out["spill_reason"] = f"too few isolated receivers ({len(I)})"
        return out
    S = np.flatnonzero(s_mask)
    lin = ex.lin
    mS = np.asarray(lin[S].mean(axis=0)).ravel()
    mI = np.asarray(lin[I].mean(axis=0)).ravel()
    ratio = mS / np.maximum(mI, 1e-3)
    cand = np.flatnonzero((ratio >= cfg.spill_marker_ratio) & (mS >= cfg.spill_marker_min_cp10k))
    excl = {ex.gi.get(ligand, -1), ex.gi.get(receptor, -1)}
    cand = np.array([c for c in cand if c not in excl], dtype=np.int64)
    if len(cand) < cfg.spill_min_markers:
        out["spill_reason"] = f"fewer than {cfg.spill_min_markers} sender markers"
        out["spill_n_markers"] = int(len(cand))
        return out
    cand = cand[np.argsort(-ratio[cand])][: cfg.spill_n_markers]
    TI = np.concatenate([T, I])
    is_t = np.zeros(len(TI), bool)
    is_t[: len(T)] = True
    Mk = np.asarray(lin[TI][:, cand].todense())
    den = mS[cand]

    def phi_of(mask_t: np.ndarray) -> np.ndarray:
        # mask_t: (b, n_TI) boolean
        nt = mask_t.sum(1, keepdims=True)
        ni = (~mask_t).sum(1, keepdims=True)
        mt = (mask_t.astype(float) @ Mk) / nt
        mi = ((~mask_t).astype(float) @ Mk) / ni
        return np.median((mt - mi) / np.maximum(den[None, :] - mi, 1e-6), axis=1)

    phi = float(phi_of(is_t[None, :])[0])
    nulls = []
    done = 0
    while done < cfg.n_perms:
        b = min(cfg.perm_chunk, cfg.n_perms - done)
        keys = np.argsort(rng.random((b, len(TI))), axis=1)
        nulls.append(phi_of(is_t[keys]))
        done += b
    null = np.concatenate(nulls)
    p_phi = _addone_p(null, phi)
    detected = bool(phi >= cfg.spill_phi_min and p_phi < cfg.alpha)
    r_lin = ex.gene(receptor)[0]
    dR_obs = float(r_lin[T].mean() - r_lin[I].mean())
    dR_exp = float(phi * (r_lin[S].mean() - r_lin[I].mean()))
    if detected and dR_obs > 0 and dR_exp >= cfg.spill_explain_frac * dR_obs:
        status = "FAIL"
    else:
        status = "PASS"
    out.update({
        "spill_status": status, "spill_phi": phi, "spill_p_phi": p_phi, "spill_detected": detected,
        "spill_n_markers": int(len(cand)), "spill_markers": ";".join(ex.genes[c] for c in cand),
        "spill_n_touching": int(len(T)), "spill_n_isolated": int(len(I)),
        "spill_dR_obs": dR_obs, "spill_dR_exp": dR_exp,
        "spill_explained_frac": float(dR_exp / dR_obs) if dR_obs > 0 else float("nan"),
        "spill_reason": ("contamination detected and explains >= {:.0%} of the touching-receiver receptor excess".format(
            cfg.spill_explain_frac) if status == "FAIL" else
            ("no contamination detected" if not detected else
             ("no receptor excess in touching receivers" if dR_obs <= 0 else
              "contamination explains < {:.0%} of the receptor excess".format(cfg.spill_explain_frac)))),
    })
    return out


def contact_statistics(ctx: SpatialContext, ex: ExpressionData, sender: str, receiver: str, ligand: str,
                       receptor: str, cfg: AuditConfig | None = None, seed: int | None = None,
                       run_spillover: bool = True) -> dict[str, Any]:
    """Compute every statistic needed by the v2 claim audit for one claim."""
    cfg = cfg or ctx.cfg
    rng = np.random.default_rng(cfg.seed if seed is None else seed)
    s_mask = ctx.labels == sender
    r_mask = ctx.labels == receiver
    L_lin, L_log, L_pos = ex.gene(ligand)
    R_lin, R_log, _ = ex.gene(receptor)
    S = np.flatnonzero(s_mask)
    R = np.flatnonzero(r_mask)
    m: dict[str, Any] = {"sender": sender, "receiver": receiver, "ligand": ligand, "receptor": receptor,
                         "n_senders": int(len(S)), "n_receivers": int(len(R)), "n_perms": cfg.n_perms,
                         "min_attainable_p": cfg.min_p, "cell_diameter_um": ctx.diameter}
    si, rj = ctx.directed_pairs(s_mask, r_mask)
    m["n_edges"] = int(len(si))
    if len(S) == 0 or len(R) == 0 or len(si) == 0:
        m["status"] = "no_edges"
        return m
    pos_s = np.full(ctx.n, -1)
    pos_s[S] = np.arange(len(S))
    pos_r = np.full(ctx.n, -1)
    pos_r[R] = np.arange(len(R))
    es, er = pos_s[si], pos_r[rj]
    Ls, Rr = L_log[S], R_log[R]
    E = float(np.mean(Ls[es] * Rr[er]))
    M = float(Ls.mean() * Rr.mean())
    Q = E / M if M > 0 else float("nan")
    m.update({"E": E, "M": M, "Q": Q})
    zs = np.zeros(len(S), np.int64)
    zr = np.zeros(len(R), np.int64)
    nullE_g = contact_nulls(Ls, Rr, es, er, zs, zr, cfg.n_perms, rng, cfg.perm_chunk)
    nullE_b = contact_nulls(Ls, Rr, es, er, ctx.tiles[S], ctx.tiles[R], cfg.n_perms, rng, cfg.perm_chunk)
    Qg, Qb = nullE_g / M, nullE_b / M
    m["p_contact_global"] = _addone_p(Qg, Q)
    m["p_contact_block"] = _addone_p(Qb, Q)
    m["z_contact_global"] = float((Q - Qg.mean()) / Qg.std(ddof=1)) if Qg.std() > 0 else float("nan")
    m["z_contact_block"] = float((Q - Qb.mean()) / Qb.std(ddof=1)) if Qb.std() > 0 else float("nan")
    # naive abundance score
    nullM = label_matched_null(L_log, R_log, len(S), len(R), cfg.n_perms, rng)
    m["p_M_label"] = _addone_p(nullM, M)
    # ring control
    exposed = np.zeros(ctx.n, bool)
    exposed[rj[L_pos[si]]] = True
    Lplus = S[L_pos[S]]
    m["n_exposed"] = int(exposed[R].sum())
    # within-tile receiver contrast (amendment 1): exposed vs unexposed receivers, Rg permuted within tiles
    xr = exposed[R]
    blk_ok = xr.sum() >= cfg.min_ring_group and (~xr).sum() >= cfg.min_ring_group
    m["block_contrast_run"] = bool(blk_ok)
    if blk_ok:
        Rvals = R_log[R]
        B_obs = float(Rvals[xr].mean() - Rvals[~xr].mean())
        nullB = ring_null(Rvals, xr, ctx.tiles[R], cfg.n_perms, rng, cfg.perm_chunk)
        m["block_contrast_diff"] = B_obs
        m["p_block_contrast"] = _addone_p(nullB, B_obs)
    ring_ok = False
    if len(Lplus) > 0:
        d_lp = KDTree(ctx.xy[Lplus]).query(ctx.xy[R], k=1)[0].ravel()
        ctrl = (~exposed[R]) & (d_lp <= cfg.ring_factor * ctx.diameter)
        xmask = exposed[R]
        m["n_ring_controls"] = int(ctrl.sum())
        if xmask.sum() >= cfg.min_ring_group and ctrl.sum() >= cfg.min_ring_group:
            ring_ok = True
            sel = xmask | ctrl
            vals = R_log[R][sel]
            is_x = xmask[sel]
            grp = ctx.tiles[R][sel]
            T_obs = float(vals[is_x].mean() - vals[~is_x].mean())
            nullT = ring_null(vals, is_x, grp, cfg.n_perms, rng, cfg.perm_chunk)
            m["ring_diff"] = T_obs
            m["p_ring"] = _addone_p(nullT, T_obs)
    else:
        m["n_ring_controls"] = 0
    m["ring_run"] = ring_ok
    # plant-and-shift power at cfg.power_effect
    R_lin_p = R_lin.copy()
    R_lin_p[exposed] *= cfg.power_effect
    R_log_p = np.log1p(R_lin_p)
    Rr_p = R_log_p[R]
    E_p = float(np.mean(Ls[es] * Rr_p[er]))
    M_p = float(Ls.mean() * Rr_p.mean())
    shift_q = E_p / M_p - Q
    pw = {"power_global": _shift_power(Qg, shift_q, cfg.alpha)}
    m["power_block_q_info"] = _shift_power(Qb, shift_q, cfg.alpha)
    if blk_ok:
        Rp = R_log_p[R]
        shift_b = float(Rp[xr].mean() - Rp[~xr].mean()) - B_obs
        pw["power_block_contrast"] = _shift_power(nullB, shift_b, cfg.alpha)
    if ring_ok:
        vals_p = R_log_p[R][sel]
        shift_t = float(vals_p[is_x].mean() - vals_p[~is_x].mean()) - T_obs
        pw["power_ring"] = _shift_power(nullT, shift_t, cfg.alpha)
    m.update(pw)
    m["power_chain"] = float(min(pw.values()))
    if run_spillover:
        m.update(spillover_check(ctx, ex, s_mask, r_mask, ligand, receptor, cfg, rng))
    m["status"] = "evaluated"
    return m
