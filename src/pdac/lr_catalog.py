"""Curated ligand-receptor catalogue for the PDAC audit.

Declared on 2026-10-24 BEFORE any PDAC data were loaded (see docs/PDAC_PROGRESS_LOG.md).
  * BREAST_CURATED_20: the 20 pairs of src/pipeline.py::LR_PAIRS (paper panel).
  * BREAST_MULTI_AXIS: the two extra breast axes of src/run_multi_axis_expansion.py.
  * PDAC_ADDITIONS: checkpoint / chemokine / stromal pairs added for the PDAC focus
    (immune-excluded stroma, CXCL12-CXCR4, checkpoint axes).
Only pairs with BOTH genes on the panel are tested; the coverage table reports all of them.
"""
from pipeline import LR_PAIRS

BREAST_CURATED_20 = list(LR_PAIRS)
BREAST_MULTI_AXIS = [("CD86", "CTLA4"), ("PTN", "SDC4")]
PDAC_ADDITIONS = [
    ("CD80", "CTLA4"), ("CD80", "CD28"), ("CD86", "CD28"),
    ("PDCD1LG2", "PDCD1"), ("LGALS9", "HAVCR2"),
    ("SPP1", "CD44"), ("FN1", "ITGB1"),
    ("CXCL13", "CXCR5"), ("CCL19", "CCR7"), ("CCL21", "CCR7"),
    ("CXCL12", "ACKR3"),
]


def catalogue() -> list[tuple[str, str, str]]:
    out = []
    seen = set()
    for src, pairs in (("breast_curated_20", BREAST_CURATED_20), ("breast_multi_axis", BREAST_MULTI_AXIS),
                       ("pdac_addition", PDAC_ADDITIONS)):
        for l, r in pairs:
            if (l, r) not in seen:
                seen.add((l, r))
                out.append((l, r, src))
    return out
