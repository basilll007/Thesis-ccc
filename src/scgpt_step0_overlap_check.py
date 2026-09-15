"""Step 0: environment check and gene-panel/vocab overlap for scGPT zero-shot baseline.

Standalone script, run inside the `scgpt-zeroshot` conda env. Does not modify
any existing results/*.csv, data/processed/*.h5ad, or src/*.py files.
"""
import json
import time
import sys

import anndata as ad
import scgpt

CKPT_DIR = r"F:\Thesis\models\scgpt_whole_human"
H5AD_PATH = r"F:\Thesis\data\processed\xenium_breast_baseline.h5ad"


def main():
    t0 = time.time()

    print("scgpt version:", getattr(scgpt, "__version__", "unknown"))
    with open(rf"{CKPT_DIR}\args.json") as f:
        args = json.load(f)
    print("checkpoint args (subset):", {k: args[k] for k in [
        "embsize", "nlayers", "nheads", "d_hid", "fast_transformer", "n_bins", "max_seq_len"
    ]})

    with open(rf"{CKPT_DIR}\vocab.json") as f:
        vocab = json.load(f)
    vocab_genes = set(vocab.keys())
    print("scGPT gene vocab size:", len(vocab_genes))

    adata = ad.read_h5ad(H5AD_PATH)
    panel_genes = list(adata.var_names)
    n_panel = len(panel_genes)
    print("n_genes_in_panel:", n_panel)

    # direct overlap
    direct_overlap = [g for g in panel_genes if g in vocab_genes]

    # case-insensitive alias check for genes missing from direct overlap
    vocab_upper_map = {}
    for vg in vocab_genes:
        vocab_upper_map.setdefault(vg.upper(), vg)

    missing_direct = [g for g in panel_genes if g not in vocab_genes]
    recovered_by_case = []
    still_missing = []
    for g in missing_direct:
        if g.upper() in vocab_upper_map:
            recovered_by_case.append((g, vocab_upper_map[g.upper()]))
        else:
            still_missing.append(g)

    overlap_count = len(direct_overlap) + len(recovered_by_case)
    overlap_fraction = overlap_count / n_panel

    print("n_genes_in_scgpt_vocab (total vocab):", len(vocab_genes))
    print("overlap_count (direct):", len(direct_overlap))
    print("overlap_count (direct + case-recovered):", overlap_count)
    print("overlap_fraction:", round(overlap_fraction, 4))
    print()
    print("genes recovered only via case-insensitive match (n={}):".format(len(recovered_by_case)))
    for g, vg in recovered_by_case:
        print(f"  panel='{g}' -> vocab='{vg}'")
    print()
    print("panel genes NOT in scGPT vocab even after case check (n={}):".format(len(still_missing)))
    for g in still_missing:
        print(f"  {g}")

    elapsed = time.time() - t0
    print()
    print(f"Step 0 runtime: {elapsed:.2f} s")

    gate_pass = overlap_fraction >= 0.5
    print("GATE (overlap_fraction >= 0.5):", "PASS" if gate_pass else "FAIL")

    # Save the final gene list actually usable (direct matches; case-aliased genes use panel name as key)
    with open(r"F:\Thesis\models\scgpt_whole_human\_step0_overlap_genes.json", "w") as f:
        json.dump({
            "n_genes_in_panel": n_panel,
            "n_genes_in_scgpt_vocab": len(vocab_genes),
            "overlap_count": overlap_count,
            "overlap_fraction": overlap_fraction,
            "direct_overlap_genes": direct_overlap,
            "case_recovered_genes": recovered_by_case,
            "still_missing_genes": still_missing,
            "gate_pass": gate_pass,
            "runtime_sec": elapsed,
        }, f, indent=2)

    return gate_pass


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
