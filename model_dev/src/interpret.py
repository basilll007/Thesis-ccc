"""Phase 5: interpretation, gated on the graph model beating the strongest baseline.

What (plan: model_dev/plans/phase5_interpretation.md): tumour-level permutation importance by neighbour class
(within-section, within-class permutation of neighbour sources) and by neighbour gene group (zeroing the group's
standardised expression in all neighbours), using the Phase 3 checkpoints on their own held-out folds. The in-silico
perturbation runs only if the cross-type contribution CI excludes 0, and its rows are labelled "model-based
hypothesis". Attention weights are never used as explanations.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mdcommon as C  # noqa: E402
import assemble as A  # noqa: E402
import features as F  # noqa: E402
import models as M  # noqa: E402
import train as T  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

PHASE = "phase5"


def class_permutation(store0, classes_arr: np.ndarray, cls_id: int, seed: int, device) -> torch.Tensor:
    """Global index map: cells of `cls_id` are permuted within their section; everything else maps to itself."""
    rng = np.random.default_rng([C.SEED, 5, seed, cls_id])
    perm = np.arange(store0.n)
    sec = store0.section_of
    for s in np.unique(sec):
        idx = np.where((sec == s) & (classes_arr == cls_id))[0]
        perm[idx] = idx[rng.permutation(len(idx))]
    return torch.as_tensor(perm).long().to(device)


def predict(model, store, g, cov, rc, cfg, device, remap=None, x_transform=None) -> np.ndarray:
    model.eval()
    out = []
    with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=device == "cuda"):
        for i in range(0, len(g), 1024):
            tb = M.TreeBatch(store, g[i:i + 1024], cov[i:i + 1024], rc, "full", cfg["k1"], cfg["k2"], remap, x_transform)
            out.append(model(tb).float().cpu().numpy())
    return np.vstack(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(C.MODEL_DEV / "configs" / "interpret.yaml"))
    args = ap.parse_args()
    cfg = C.load_config(args.config)
    log = C.start_run(PHASE, "interpret", cfg)
    summ3 = pd.read_csv(C.RESULTS / "phase3" / "graph_summary.csv")
    sel = json.loads((C.RESULTS / "phase3" / "architecture_selection.json").read_text())
    arch, lr = sel["selected_arch"], float(sel["selected_lr"])
    gated = summ3[(summ3.arch == arch) & summ3.beats_baseline.astype(str).eq("True") & summ3.target.str.contains(":")]
    gate = {"selected_arch": arch, "targets_passing_gate": gated[["receiver", "target", "gain_vs_baseline", "gain_ci_low"]].to_dict("records"),
            "rule": "graph model beats strongest baseline (paired tumour-bootstrap CI of gain excludes 0)"}
    if gated.empty:
        gate["decision"] = "Phase 5 not run: no target where the graph model beats the strongest baseline"
        C.save_json(gate, PHASE, "gate.json", log)
        log.info(gate["decision"])
        return
    gate["decision"] = "Phase 5 run for the targets listed"
    C.save_json(gate, PHASE, "gate.json", log)
    abl = pd.read_csv(C.RESULTS / "phase3" / "ablation_summary.csv")
    gcfg = C.load_config(C.MODEL_DEV / "configs" / "graph.yaml")
    if not torch.cuda.is_available():
        raise RuntimeError("GPU-only policy: CUDA not available")
    device = "cuda"
    fold_table = A.folds("grouped"); n_folds = int(fold_table.fold.max()) + 1
    spec = A.targets_spec()
    programs = json.loads((C.REPO_ROOT / "results/chain/phase1b/programs_used.json").read_text())
    genes = [str(g) for g in A.load_section(A.section_files("development")[0])["genes"]]
    gene_groups = {k: [genes.index(g) for g in v if g in genes] for k, v in {**programs["programs"], **programs["lineage"]}.items()}
    gene_groups["other_G_genes"] = sorted(set(range(len(genes))) - {i for v in gene_groups.values() for i in v})
    store0, offsets, classes = T.build_store("development", max(gcfg["k1"], gcfg["k2"]), device, log)
    classes_arr = store0.cls.cpu().numpy()
    rows, pert_rows = [], []
    for receiver in sorted(gated.receiver.unique()):
        s = spec[receiver]
        data = A.build_receiver("development", receiver, log, blocks=["own", "tile"])
        tnames = data.target_names
        targets = gated[gated.receiver == receiver].target.tolist()
        rc = classes.index(s["class"])
        g_all = torch.as_tensor(np.array([offsets[x] + r for x, r in zip(data.section, data.row)])).long().to(device)
        raw = data.matrix(["own", "tile"])
        drops = {}   # (group, target) -> {tumour: [per-seed drop]}
        changes = {}
        for k in range(n_folds):
            tr_in, val, tr_all, test = A.fold_masks(data, fold_table, k, n_folds)
            te = np.where(test)[0]
            ytr_mean = data.Y[tr_all].mean(axis=0)
            for seed in gcfg["seeds"]:
                ck = torch.load(C.CACHE / "checkpoints" / "phase3" / f"graph_{receiver}_full_{arch}_lr{lr}_seed{seed}_fold{k}.pt", weights_only=False)
                store = T.StandardisedStore(store0, ck["mu"].to(device), ck["sd"].to(device), ck["cmu"].to(device), ck["csd"].to(device))
                cov = torch.as_tensor(F.apply_scaler(raw, ck["smu"], ck["ssd"])).to(device)
                model = M.TreeModel(store.X.shape[1], len(classes), raw.shape[1], len(tnames), gcfg["hidden"], gcfg["emb"],
                                    "attention" if arch == "attention" else "mean", gcfg["dropout"], gcfg["heads"]).to(device)
                model.load_state_dict(ck["state"])
                g, c = g_all[te], cov[te]
                base = predict(model, store, g, c, rc, gcfg, device)
                groups = {f"class:{cn}": ("remap", class_permutation(store0, classes_arr, ci, seed, device)) for ci, cn in enumerate(classes)}
                for gn, gi in gene_groups.items():
                    gm = torch.ones(len(genes), device=device); gm[torch.as_tensor(gi).to(device)] = 0
                    groups[f"genes:{gn}"] = ("xt", lambda x, cls, gm=gm: x * gm)
                for gname, (kind, obj) in groups.items():
                    pred = predict(model, store, g, c, rc, gcfg, device, remap=obj if kind == "remap" else None,
                                   x_transform=obj if kind == "xt" else None)
                    for t in targets:
                        j = tnames.index(t)
                        for tum in np.unique(data.group[te]):
                            m = data.group[te] == tum
                            sse_b = float(((data.Y[te][m, j] - base[m, j]) ** 2).sum())
                            sse_p = float(((data.Y[te][m, j] - pred[m, j]) ** 2).sum())
                            sst = float(((data.Y[te][m, j] - ytr_mean[j]) ** 2).sum())
                            drops.setdefault((gname, t), {}).setdefault(tum, []).append((sse_p - sse_b, sst))
                # in-silico perturbation (double gate): cross-type contribution CI excludes 0 for this target
                for t in targets:
                    a = abl[(abl.receiver == receiver) & (abl.target == t) & (abl["mode"] == "same_type")]
                    if a.empty or not bool(a.diff_ci_excludes_0.iloc[0]):
                        continue
                    for pname, gi in gene_groups.items():
                        if ":" not in pname or pname.split(":")[0] == {v: k for k, v in F.TYPED.items()}.get(s["class"]):
                            continue  # only other-type programs
                        ptype = F.TYPED[pname.split(":")[0]]
                        ci = classes.index(ptype)
                        for q in cfg["quantiles"]:
                            qv = torch.as_tensor(np.quantile(store0.X[store0.cls == ci].float().cpu().numpy()[:, gi], q, axis=0)).to(device)
                            qv = (qv - ck["mu"][gi].to(device)) / ck["sd"][gi].to(device)
                            gi_t = torch.as_tensor(gi).to(device)

                            def xt(x, cls, qv=qv, gi_t=gi_t, ci=ci):
                                x = x.clone()
                                sel_c = (cls == ci)
                                x[..., gi_t] = torch.where(sel_c.unsqueeze(-1), qv.expand(x[..., gi_t].shape), x[..., gi_t])
                                return x
                            pred = predict(model, store, g, c, rc, gcfg, device, x_transform=xt)
                            j = tnames.index(t)
                            for tum in np.unique(data.group[te]):
                                m = data.group[te] == tum
                                changes.setdefault((pname, q, t), {}).setdefault(tum, []).append(float((pred[m, j] - base[m, j]).mean()))
                del model
            log.info("%s fold %d importance done", receiver, k)
        for (gname, t), per_t in drops.items():
            tums = sorted(per_t)
            d = np.array([np.mean([x[0] for x in per_t[u]]) for u in tums]); sst = np.array([per_t[u][0][1] for u in tums])
            point = d.sum() / sst.sum()
            rng = np.random.default_rng(C.SEED); idx = rng.integers(0, len(tums), size=(cfg["n_boot"], len(tums)))
            boots = d[idx].sum(1) / sst[idx].sum(1)
            rows.append({"receiver": receiver, "target": t, "group": gname, "r2_drop": point, "ci_low": float(np.quantile(boots, .025)),
                         "ci_high": float(np.quantile(boots, .975)), "n_tumours": len(tums)})
        for (pname, q, t), per_t in changes.items():
            tums = sorted(per_t); v = np.array([np.mean(per_t[u]) for u in tums])
            rng = np.random.default_rng(C.SEED); idx = rng.integers(0, len(tums), size=(cfg["n_boot"], len(tums)))
            boots = v[idx].mean(1)
            pert_rows.append({"receiver": receiver, "target": t, "neighbour_program": pname, "quantile": q,
                              "predicted_change": float(v.mean()), "ci_low": float(np.quantile(boots, .025)),
                              "ci_high": float(np.quantile(boots, .975)), "label": "model-based hypothesis"})
    imp = pd.DataFrame(rows).sort_values(["receiver", "target", "r2_drop"], ascending=[True, True, False])
    C.save_table(imp, PHASE, "permutation_importance.csv", log)
    C.save_table(pd.DataFrame(pert_rows), PHASE, "in_silico_perturbation.csv", log)
    log.info("top importances:\n%s", imp.groupby(["receiver", "target"]).head(8).round(4).to_string(index=False))
    log.info("done")


if __name__ == "__main__":
    main()
