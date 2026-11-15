# v3 Phase 4 plan — replication on GSE310352 (GATED)

Written 2026-11-15. Gate: the file `LICENSE_CONFIRMED_GSE310352.txt` at the repository root must exist
(created by the user after checking the data licence). Checked at run time by `src/pdac/v3_phase4_gate.py`.

If the gate opens: data = `data/pdac/geo_gse310352/` (loader `pdac_loaders.load_geo_gse310352`, 8 slides; coordinates in px,
um/px UNVERIFIED -> block tiles would be defined in px with the Phase 1 sizes converted only if a verified um/px exists;
otherwise block nulls are run in units of the median nearest-neighbour distance (5x, 10x, 20x) and labelled as such).
Labels: marker scoring and Leiden consensus (Phase 3 `annotate`), no author labels. Nulls: Phase 1 set. Report per axis:
same direction as Mendeley (yes/no), background flags from negative probes, and whether the 8 slides are from patients
different from the Mendeley set (needs the GEO sample metadata; otherwise UNVERIFIED and the word "replication" is not used).

If the gate is closed: log BLOCKED (licence unconfirmed) and skip. Status at plan time: file absent -> BLOCKED.
