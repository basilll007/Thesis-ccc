# v2 Phase 6 plan: cellalibi verdicts on all real sections

Committed before any Phase 6 audit is run. The engine is the frozen Phase 4 engine (`cellalibi` at fc9bab8 for
`contact.py` / `verdicts.py`). The only additions are loaders and program features; they do not change engine
behaviour.

## Data

- **Mendeley PDAC-SMI**: all 16 sections, thesis `consensus_T` labels (the v3 primary annotation).
  - Loader: `cellalibi.datasets.load_mendeley_section`.
  - Coordinates in um. Tumour (patient) = thesis section-to-tumour map, giving 13 tumours.
- **GSE310352**: all 8 slides, thesis `marker` labels (the only annotation with all four focus types on 8/8 slides,
  per thesis claim 17).
  - Loader: `cellalibi.datasets.load_gse310352_slide`, using the thesis chain cache for cells and labels.
  - Coordinates in px (um/px UNVERIFIED), so the tile size is set to 10 x d (median nearest-neighbour distance). On
    Mendeley the default 100 um tile is 10.2 d. The ring radius (4 d) and spillover distances are already relative
    to d.
  - Patient independence of slides is UNVERIFIED (thesis claim 13), so each slide is its own unit and GSE310352 is
    never pooled with Mendeley.

## Claims

1. **Catalogue L-R claims**:
   - Every thesis catalogue pair (`src/pdac/lr_catalog.py`, 33 pairs) with both genes on the panel;
   - for every ordered pair of distinct focus types among epithelial_ductal_tumor, CAF_fibroblast, T_cell and
     myeloid (12 ordered pairs).
   - Claims whose sender or receiver type is absent are recorded as NOT_TESTABLE (no edges).
2. **Receiver-program claims**: "contact with sender type S raises receiver program P in receiver type R".
   - Programs use the thesis chain gene lists (`src/chain/chain_phase1_cell_states.py::PROGRAMS`).
   - Program score = mean CP10K of the program genes on the panel (linear, so spillover mixing algebra holds).
   - Ligand = "any S contact" (every S cell counts as positive).
   - Claims:
     - CAF_fibroblast -> myeloid: C1Q_TAM, SPP1_TAM;
     - CAF_fibroblast -> epithelial_ductal_tumor: classical, basal_like;
     - myeloid -> epithelial_ductal_tumor: classical, basal_like;
     - T_cell -> epithelial_ductal_tumor: classical, basal_like.

Engine settings: package defaults (1000 permutations, alpha 0.05, power floor 0.5 at 1.5x, etc.). No multiplicity
correction inside the engine. The false-SUPPORTED rate per claim under the null comes from Phase 4 S0 and is quoted
next to the count of SUPPORTED claims.

## Aggregation (patient is the unit of inference)

- **Per tumour**: the primary verdict across that tumour's sections is the most frequent one. Ties go to the
  verdict that comes first in the priority order (NOT_TESTABLE > SPILLOVER > ABUNDANCE/NOT_SUPPORTED > ZONATION >
  SUPPORTED), which is the conservative choice.
- **Per claim, cohort level**:
  - SUPPORTED if SUPPORTED in a majority of the tumours where it is testable, and testable in >= 3 tumours;
  - NOT_TESTABLE if testable in < 3 tumours;
  - otherwise the most frequent non-SUPPORTED verdict among testable tumours.
  - Mendeley and GSE310352 are reported separately; slides are the units for GSE310352.
- Artefact frequency: how often each check FAILs, and how often each primary verdict occurs, over all
  claim x section audits.

## Comparison with recorded thesis conclusions

Each claim the thesis reports drew a conclusion about is compared:
- v3 claims register items 1-17;
- model-v2 CAF->macrophage C1Q_TAM = EXPLAINED_BY_SPILLOVER;
- the v1 demo verdicts;
- FINDINGS.md entries.

Every mismatch is explained. If an explanation cannot be found, it is marked "unexplained".

## Runtime and memory

Per section: wall time and peak RSS of the worker for load + all claims, via psutil. Free RAM is logged before start.
- 4 workers, 1 BLAS thread each.
- One resumable CSV chunk per section.
