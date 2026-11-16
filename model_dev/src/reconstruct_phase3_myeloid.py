"""One-off reconstruction of Phase 3 myeloid bookkeeping after the tumour OOM (logged in PROGRESS_LOG)."""
import sys, re, glob
sys.path.insert(0, 'model_dev/src')
import numpy as np, pandas as pd, mdcommon as C, assemble as A
ph = C.RESULTS / 'phase3'
folds = A.folds('grouped'); t2f = dict(zip(folds.tumour, folds.fold))
per = pd.read_csv(ph / 'per_tumour_sse_graph_myeloid.csv'); per['fold'] = per.tumour.astype(str).map(t2f)
spec = A.targets_spec()['myeloid']; nprog = len(spec['programs'])
names = spec['programs'] + spec['genes']
rows, comp = [], []
log = open(ph / 'run_train_20261116T012817.log', encoding='utf-8').read()
pat = re.compile(r"myeloid full (\w+) lr=([\d.e-]+) seed=(\d) fold (\d): (\d+) epochs, val mse ([\d.]+), (\d+)s, peak (\d+) MiB")
for m in pat.finditer(log):
    arch, lr, seed, k, ep, vm, wall, peak = m.groups(); lr = float(lr); seed = int(seed); k = int(k)
    cur = pd.read_csv(ph / 'curves' / f'graph_myeloid_full_{arch}_lr{lr}_seed{seed}_fold{k}.csv')
    assert len(cur) == int(ep) and abs(cur.val_mse.min() - float(vm)) < 1e-3
    d = per[(per.arch == arch) & (per.lr == lr) & (per.seed == seed) & (per.fold == k)]
    for t, x in d.groupby('target'):
        kind = 'gene_vector' if t == 'genes_pooled' else ('program' if t in spec['programs'] else 'gene')
        rows.append({'receiver': 'myeloid', 'mode': 'full', 'arch': arch, 'lr': lr, 'seed': seed, 'fold': k, 'target': t,
                     'target_kind': kind, 'inner_val_mse': float(cur.val_mse.min()), 'epochs': int(ep),
                     'test_r2_fold': 1 - x.sse.sum() / x.sst.sum()})
    comp.append({'receiver': 'myeloid', 'mode': 'full', 'arch': arch, 'lr': lr, 'seed': seed, 'fold': k, 'wall_s': float(wall),
                 'epochs': int(ep), 'peak_vram_mib': float(peak), 'device': 'cuda', 'n_train': np.nan,
                 'checkpoint': f'model_cache/checkpoints/phase3/graph_myeloid_full_{arch}_lr{lr}_seed{seed}_fold{k}.pt'})
print(len(rows), len(comp))
assert len(comp) == 60
pd.DataFrame(rows).to_csv(ph / 'results_graph_myeloid.csv', index=False)
pd.DataFrame(comp).to_csv(ph / 'compute_graph_myeloid.csv', index=False)
