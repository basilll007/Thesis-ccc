"""Step 1: zero-shot scGPT cell embeddings vs raw-expression PCA baseline gate.

Standalone script for the `scgpt-zeroshot` conda env. Reads
data/processed/xenium_breast_baseline.h5ad (read-only) and writes only to
results/embedding_baseline/. Does not modify any existing results/src files.
"""
import json
import os
import time

import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA

# Windows has no os.sched_getaffinity; scgpt's DataLoader worker-count logic needs it.
# Return an empty set so scgpt's num_workers=min(len(affinity), batch_size) resolves to 0,
# which avoids Windows spawn-based multiprocessing (the DataLoader's Dataset is a local
# class defined inside get_batch_cell_embeddings and cannot be pickled to a worker process).
if not hasattr(os, "sched_getaffinity"):
    os.sched_getaffinity = lambda pid: set()

import scgpt.tasks.cell_emb as ce

CKPT_DIR = r"F:\Thesis\models\scgpt_whole_human"
H5AD_PATH = r"F:\Thesis\data\processed\xenium_breast_baseline.h5ad"
OUT_DIR = r"F:\Thesis\results\embedding_baseline"
SUBSET_TEST = False  # set True to smoke-test on a small subset first


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    t0 = time.time()

    adata = ad.read_h5ad(H5AD_PATH)
    print("Loaded adata:", adata.shape)

    if SUBSET_TEST:
        adata = adata[:50].copy()
        print("SUBSET_TEST active, using", adata.shape)

    adata_emb_input = adata.copy()

    cell_emb_adata = ce.embed_data(
        adata_emb_input,
        model_dir=CKPT_DIR,
        gene_col="index",
        max_length=1200,
        batch_size=64,
        obs_to_save=["cell_type"],
        device="cpu",
        use_fast_transformer=False,
        return_new_adata=False,
    )

    scgpt_emb = np.asarray(cell_emb_adata.obsm["X_scGPT"])
    print("scGPT embedding shape:", scgpt_emb.shape)
    n_genes_used = cell_emb_adata.n_vars
    print("n_genes actually used by scGPT (post vocab filter):", n_genes_used)

    t1 = time.time()
    step1_embed_runtime = t1 - t0
    print(f"scGPT zero-shot embedding extraction runtime: {step1_embed_runtime:.2f} s")

    # --- PCA baseline on the SAME gene-restricted expression matrix ---
    # Restrict to the exact same genes scGPT actually used (post vocab overlap filter).
    used_genes = list(cell_emb_adata.var_names)
    adata_pca_input = adata[:, used_genes].copy()
    X = adata_pca_input.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    n_pcs = min(50, X.shape[1] - 1, X.shape[0] - 1)
    pca = PCA(n_components=n_pcs, random_state=42)
    pca_emb = pca.fit_transform(X)
    print("PCA baseline embedding shape:", pca_emb.shape, "(n_pcs =", n_pcs, ")")

    # --- silhouette scores against existing cell_type labels ---
    labels = adata.obs["cell_type"].astype(str).values

    sil_scgpt = silhouette_score(scgpt_emb, labels)
    sil_pca = silhouette_score(pca_emb, labels)

    print(f"Silhouette (scGPT zero-shot embeddings): {sil_scgpt:.4f}")
    print(f"Silhouette (PCA baseline, {n_pcs} PCs, same {n_genes_used} genes): {sil_pca:.4f}")

    # Gate: scGPT must be within ~10% of PCA (allow it to also exceed PCA)
    if sil_pca > 0:
        ratio = sil_scgpt / sil_pca
    else:
        ratio = float("inf") if sil_scgpt >= 0 else float("-inf")
    gate_pass = sil_scgpt >= 0.90 * sil_pca
    print(f"ratio scGPT/PCA = {ratio:.4f}")
    print("GATE (scGPT silhouette within ~10% of, or better than, PCA):", "PASS" if gate_pass else "FAIL")

    t2 = time.time()

    # --- UMAP figures ---
    def make_umap_adata(emb, name):
        a = sc.AnnData(X=emb)
        a.obs["cell_type"] = adata.obs["cell_type"].values
        sc.pp.neighbors(a, use_rep="X", random_state=42)
        sc.tl.umap(a, random_state=42)
        return a

    umap_scgpt = make_umap_adata(scgpt_emb, "scgpt")
    umap_pca = make_umap_adata(pca_emb, "pca")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    cell_types = sorted(adata.obs["cell_type"].astype(str).unique())
    palette = dict(zip(cell_types, plt.cm.tab10.colors))

    for ax, umap_a, title, sil in [
        (axes[0], umap_scgpt, "scGPT zero-shot embeddings", sil_scgpt),
        (axes[1], umap_pca, f"PCA baseline ({n_pcs} PCs, same genes)", sil_pca),
    ]:
        for ct in cell_types:
            mask = umap_a.obs["cell_type"].values == ct
            ax.scatter(
                umap_a.obsm["X_umap"][mask, 0],
                umap_a.obsm["X_umap"][mask, 1],
                s=4, alpha=0.6, label=ct, color=palette[ct],
            )
        ax.set_title(f"{title}\nsilhouette (cell_type) = {sil:.4f}")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
    axes[0].legend(markerscale=3, fontsize=8, loc="best")

    plt.tight_layout()
    fig_path_png = os.path.join(OUT_DIR, "umap_scgpt_vs_pca_baseline.png")
    fig_path_pdf = os.path.join(OUT_DIR, "umap_scgpt_vs_pca_baseline.pdf")
    plt.savefig(fig_path_png, dpi=150)
    plt.savefig(fig_path_pdf)
    plt.close(fig)
    print("Saved:", fig_path_png, fig_path_pdf)

    t3 = time.time()
    step1_total_runtime = t3 - t0

    result = {
        "n_cells": int(adata.n_obs),
        "n_genes_used_by_scgpt": int(n_genes_used),
        "n_pcs_pca_baseline": int(n_pcs),
        "silhouette_scgpt": float(sil_scgpt),
        "silhouette_pca_baseline": float(sil_pca),
        "ratio_scgpt_over_pca": float(ratio) if np.isfinite(ratio) else None,
        "gate_pass": bool(gate_pass),
        "runtime_sec_embedding_extraction": float(step1_embed_runtime),
        "runtime_sec_total_step1": float(step1_total_runtime),
        "checkpoint": "scGPT whole-human (bowang-lab)",
        "use_fast_transformer": False,
        "device": "cpu",
    }
    with open(os.path.join(OUT_DIR, "step1_silhouette_results.json"), "w") as f:
        json.dump(result, f, indent=2)

    # Also persist the scGPT embeddings for Step 2 (if gate passes) without touching original h5ad.
    np.save(os.path.join(OUT_DIR, "scgpt_cell_embeddings.npy"), scgpt_emb)
    with open(os.path.join(OUT_DIR, "scgpt_used_genes.json"), "w") as f:
        json.dump(used_genes, f)

    print(json.dumps(result, indent=2))
    print(f"Step 1 total runtime: {step1_total_runtime:.2f} s")

    return gate_pass


if __name__ == "__main__":
    ok = main()
    print("STEP1_GATE_PASS=", ok)
