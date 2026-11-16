"""Chain-test engine (Phase 3 plan): exposures, hop regressions, mediation, block nulls with random tile offsets.

Shared by Phase 3 (real data), Phase 4 (simulation) and Phase 6 (GSE310352).
All regressions are OLS on within-section z-scored variables (standardised coefficients).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chain_features as cf  # noqa: E402

import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402

STATS = ["beta_CB", "beta_BA", "beta_total", "beta_direct", "b_mediator", "a_path", "indirect",
         "n_hop1", "n_hop2", "n_chain"]
MIN_ROWS = 30
CLASS_TO_KEY = {v: k for k, v in cf.TYPES.items()}


def within_group_perm_index(group: np.ndarray, rng) -> np.ndarray:
    """Index array p such that x[p] permutes x independently within each group (v3 lexsort construction)."""
    pos = np.argsort(group, kind="stable")
    order = np.lexsort((rng.random(group.shape[0]), group))
    p = np.empty(group.shape[0], np.int64)
    p[pos] = order
    return p


class SectionCtx:
    """Fixed per-section quantities (geometry and per-cell covariates)."""

    def __init__(self, xy, W: sp.csr_matrix, R: sp.csr_matrix, total_counts, lineage, lineage_names, tile=cf.TILE):
        self.xy = xy
        self.W = W.tocsr().astype(np.float32)
        self.n = xy.shape[0]
        self.logtot = np.log1p(np.asarray(total_counts, float))
        self.dens = np.asarray(R.sum(axis=1)).ravel().astype(float)
        self.tiles = cf.tile_codes(xy, tile)
        self.tile_inv = np.unique(self.tiles, return_inverse=True)[1]
        self.n_tiles = int(self.tile_inv.max()) + 1
        self.tile_size = tile
        self.lineage = {str(n).split(":")[1]: lineage[:, j].astype(float) for j, n in enumerate(lineage_names)}

    # ---- nulls
    def expr_perm(self, lab, rng):
        off = rng.uniform(0, self.tile_size, size=2)
        t = cf.tile_codes(self.xy, self.tile_size, off)
        _, li = np.unique(lab, return_inverse=True)
        return within_group_perm_index(t * 64 + li, rng)

    def label_perm(self, rng):
        off = rng.uniform(0, self.tile_size, size=2)
        return within_group_perm_index(cf.tile_codes(self.xy, self.tile_size, off), rng)


def _tile_ctrl(ctx, values, member, loo):
    """Tile (fixed 200 um grid) mean of `values` over members; leave-one-out for members when loo=True.

    Empty -> section mean of members, with an indicator. Uses the precomputed tile index (bincount)."""
    g = ctx.tile_inv
    v = np.where(member, values, 0.0)
    mf = member.astype(float)
    s = np.bincount(g, weights=v, minlength=ctx.n_tiles)[g]
    c = np.bincount(g, weights=mf, minlength=ctx.n_tiles)[g]
    if loo:
        s, c = s - v, c - mf
    has = c > 0
    fill = values[member].mean() if member.any() else 0.0
    return np.where(has, s / np.maximum(c, 1), fill), (~has).astype(float)


def _ols(y, X, k):
    """Standardised OLS; returns the first k coefficients (exposures). Constant controls are dropped."""
    if y.shape[0] < MIN_ROWS:
        return [np.nan] * k
    sd = X.std(axis=0)
    if np.any(sd[:k] == 0) or y.std() == 0:
        return [np.nan] * k
    keep = np.r_[np.arange(k), k + np.where(sd[k:] > 1e-12)[0]]
    Z = (X[:, keep] - X[:, keep].mean(axis=0)) / sd[keep]   # centred: no intercept needed
    yz = (y - y.mean()) / y.std()
    G = Z.T @ Z
    try:
        coef = np.linalg.solve(G, Z.T @ yz)
    except np.linalg.LinAlgError:
        coef = np.linalg.lstsq(Z, yz, rcond=None)[0]
    return list(coef[:k])


def chain_stats(ctx: SectionCtx, lab, S, col, A, B, C, triples, spill=False, geom=None, return_geom=False):
    """Statistics for one chain variant. lab: class labels; S: n x P scores; col: program -> column index.

    A, B, C are class names; triples: list of (prog_A, prog_B, prog_C) score columns to use (program names).
    Returns {triple: {stat: value}}.
    """
    W = ctx.W
    out = {}
    pa_s = sorted({t[0] for t in triples})
    pb_s = sorted({t[1] for t in triples})
    if geom is None:
        mA, mB, mC = lab == A, lab == B, lab == C
        if mA.sum() == 0 or mB.sum() == 0 or mC.sum() == 0:
            res = {t: {s: np.nan for s in STATS} for t in triples}
            return (res, None) if return_geom else res
        cntB = np.asarray(W @ mB.astype(np.float32)).ravel()
        cntA1 = np.asarray(W @ mA.astype(np.float32)).ravel()
        cidx = np.where(mC & (cntB > 0))[0]
        bidx = np.where(mB & (cntA1 > 0))[0]
        WC = W[cidx]
        P = (WC.multiply(mB[None, :].astype(np.float32))).tocsr() @ (W.multiply(mA[None, :].astype(np.float32))).tocsr()
        P = (P - P.multiply(WC)).tocsr()
        P.eliminate_zeros()
        cntA2 = np.asarray(P.sum(axis=1)).ravel()
        geom = {"mA": mA, "mB": mB, "mC": mC, "cidx": cidx, "bidx": bidx, "P": P, "cntA2": cntA2}
    mA, mB, mC, cidx, bidx, P, cntA2 = (geom[k] for k in ("mA", "mB", "mC", "cidx", "bidx", "P", "cntA2"))
    EB, _ = cf.neighbour_mean(W, S[:, [col[p] for p in pb_s]], mB)
    EA1, _ = cf.neighbour_mean(W, S[:, [col[p] for p in pa_s]], mA)
    with np.errstate(invalid="ignore", divide="ignore"):
        EA2 = (P @ S[:, [col[p] for p in pa_s]]) / cntA2[:, None]
    c2 = cntA2 > 0
    ctrl_cache = {}

    def tc(p, member, loo):
        key = (p, id(member), loo)
        if key not in ctrl_cache:
            ctrl_cache[key] = _tile_ctrl(ctx, S[:, col[p]], member, loo)
        return ctrl_cache[key]

    def spill_cols(outcome_cls, others, rows):
        if not spill:
            return []
        cols = []
        for o in others:
            k = CLASS_TO_KEY.get(o)
            if k in ctx.lineage:
                cols.append(ctx.lineage[k][rows])
        return cols

    hop1, hop2 = {}, {}
    for (pa, pb, pc) in triples:
        r = {}
        jb, ja = pb_s.index(pb), pa_s.index(pa)
        # hop 1: C <- B
        if (pb, pc) not in hop1:
            yc = S[cidx, col[pc]]
            tmc, emc = tc(pc, mC, True)
            tmb, emb = tc(pb, mB, False)
            X = np.column_stack([EB[cidx, jb], ctx.logtot[cidx], ctx.dens[cidx], tmc[cidx], emc[cidx], tmb[cidx],
                                 emb[cidx]] + spill_cols(C, [A, B], cidx))
            hop1[(pb, pc)] = (_ols(yc, X, 1)[0], len(cidx))
        # hop 2: B <- A
        if (pa, pb) not in hop2:
            yb = S[bidx, col[pb]]
            tmb_l, emb_l = tc(pb, mB, True)
            tma, ema = tc(pa, mA, False)
            X = np.column_stack([EA1[bidx, ja], ctx.logtot[bidx], ctx.dens[bidx], tmb_l[bidx], emb_l[bidx], tma[bidx],
                                 ema[bidx]] + spill_cols(B, [A, C], bidx))
            hop2[(pa, pb)] = (_ols(yb, X, 1)[0], len(bidx))
        r["beta_CB"], r["n_hop1"] = hop1[(pb, pc)]
        r["beta_BA"], r["n_hop2"] = hop2[(pa, pb)]
        # chain
        rows = cidx[c2]
        xa2, xb = EA2[c2, ja], EB[rows, jb]
        y = S[rows, col[pc]]
        tmc, emc = tc(pc, mC, True)
        tmb, emb = tc(pb, mB, False)
        tma, ema = tc(pa, mA, False)
        ctrl = [ctx.logtot[rows], ctx.dens[rows], tmc[rows], emc[rows], tmb[rows], emb[rows], tma[rows], ema[rows]]
        ctrl += spill_cols(C, [A, B], rows)
        tot = _ols(y, np.column_stack([xa2] + ctrl), 1)[0]
        dirc, bm = _ols(y, np.column_stack([xa2, xb] + ctrl), 2)
        a = _ols(xb, np.column_stack([xa2] + ctrl), 1)[0]
        r.update({"beta_total": tot, "beta_direct": dirc, "b_mediator": bm, "a_path": a, "indirect": a * bm,
                  "n_chain": len(rows)})
        out[(pa, pb, pc)] = r
    return (out, geom) if return_geom else out


def edge_count(W: sp.csr_matrix, lab, t1, t2) -> int:
    m1, m2 = (lab == t1).astype(np.float32), (lab == t2).astype(np.float32)
    return int(m1 @ (W @ m2))
