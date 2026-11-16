from common import *
import ast,re,gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc
from scipy import sparse
from scipy.stats import wilcoxon,PermutationMethod

PRIMARY='Primary PDAC'
SOURCE_SCRIPT=ROOT/'src/pdac/v3_phase5_programs.py'
_tree=ast.parse(SOURCE_SCRIPT.read_text(encoding='utf-8'))
_ns={'np':np,'pd':pd,'wilcoxon':wilcoxon,'PermutationMethod':PermutationMethod,'bh':bh,'PRIMARY':PRIMARY}
for node in _tree.body:
 if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='PROGRAMS' for t in node.targets):exec(compile(ast.Module(body=[node],type_ignores=[]),str(SOURCE_SCRIPT),'exec'),_ns)
 if isinstance(node,ast.FunctionDef) and node.name=='compare':exec(compile(ast.Module(body=[node],type_ignores=[]),str(SOURCE_SCRIPT),'exec'),_ns)
PROGRAMS=_ns['PROGRAMS'];HYPO=PROGRAMS['Hypoxia']
HOST={PRIMARY:['PRSS1','CPA1','PNLIP','AMY2A'],'Liver metastasis':['ALB','APOA1','APOA2','TTR'],'Lung metastasis':['SFTPA1','SFTPA2','SFTPB','SFTPC','SCGB1A1'],'Peritoneal metastasis':['UPK3B','CALB2','WT1','MSLN']}

def paired(frame,scenario):
 site=frame.groupby(['patient','site'],as_index=False).score.mean();site['arm']=np.where(site.site==PRIMARY,'primary','metastasis')
 arms=site.groupby(['patient','arm'],as_index=False).score.mean()
 wide=arms.pivot(index='patient',columns='arm',values='score').reindex(columns=['primary','metastasis']).dropna()
 diff=np.round((wide.metastasis-wide.primary).to_numpy(),12)
 if len(diff) and np.any(diff):test=wilcoxon(diff,method=PermutationMethod(n_resamples=np.inf));p=float(test.pvalue)
 elif len(diff):p=1.
 else:p=np.nan
 result={'scenario':scenario,'n_pairs':len(diff),'median_paired_effect':float(np.median(diff)) if len(diff) else np.nan,'p':p,'n_slides':frame.gsm.nunique(),'status':'TESTABLE' if len(diff)>=5 else 'INSUFFICIENT_PAIRS'}
 deltas=wide.assign(difference=diff).reset_index();deltas['scenario']=scenario
 return result,deltas

def score(ad,genes,key='hypoxia'):
 present=[g for g in genes if g in ad.var_names]
 if len(present)<5:raise ValueError('Fewer than five genes: '+key)
 sc.tl.score_genes(ad,present,score_name=key,ctrl_size=50,n_bins=25,ctrl_as_ref=True,random_state=42,use_raw=False)
 return ad.obs[key].to_numpy().copy()

def load(sample):
 ad=sc.read_10x_h5(SOURCE/'eda_gse277783/extracted'/sample.h5_file);ad.var_names_make_unique()
 counts=np.asarray(ad.X.sum(axis=1)).ravel();return ad,counts

def normalize(ad):
 sc.pp.normalize_total(ad,target_sum=1e4);sc.pp.log1p(ad);return ad

def coefficients(ad,sets):
 # Exact pinned Scanpy 1.10.3 bin/control selection, vectorised only at final multiplication.
 names=pd.Index(ad.var_names);pool=names.astype('string')
 means=pd.Series(np.asarray(ad.X.sum(axis=0,dtype=np.float64)).ravel()/ad.n_obs,index=pool)
 means=means[np.isfinite(means)];n_items=int(np.round(len(means)/24));cuts=means.rank(method='min')//n_items
 rr=[];cc=[];vv=[]
 for k,genes in enumerate(sets):
  gl=pd.Index(genes).intersection(names);rng=np.random.RandomState(42);controls=pd.Index([],dtype='string')
  for cut in np.unique(cuts.loc[gl]):
   candidates=cuts[cuts==cut].index
   chosen=candidates.to_series().sample(50,random_state=rng).index if len(candidates)>50 else candidates
   controls=controls.union(chosen.difference(gl))
  if not len(controls):raise ValueError('No matched controls')
  gi=names.get_indexer(gl);ci=names.get_indexer(controls)
  rr.extend(gi);rr.extend(ci);cc.extend([k]*(len(gi)+len(ci)));vv.extend([1/len(gi)]*len(gi));vv.extend([-1/len(ci)]*len(ci))
 return sparse.csc_matrix((vv,(rr,cc)),shape=(len(names),len(sets)))

def main():
 out,log=start('phase1');phase='phase1'
 try:
  meta=pd.read_csv(ROOT/'results/pdac/v3_phase5/verified_sample_metadata.csv');qc=pd.read_csv(ROOT/'results/pdac/v3_phase5/slide_qc.csv')
  samples=meta[meta.gsm.isin(qc.loc[qc.included,'gsm'])].copy().sort_values('gsm')
  base=[];reference=[]
  for s in samples.itertuples(index=False):
   spot=pd.read_csv(ROOT/'results/pdac/v3_phase5/slides'/f'{s.gsm}_spot_scores.csv.gz')
   for name in PROGRAMS:base.append({'gsm':s.gsm,'patient':s.patient,'site':s.site,'program':name,'score':float(spot[name].median())})
  base=pd.DataFrame(base);ps=base.groupby(['patient','site','program'],as_index=False).score.mean()
  _,rep=_ns['compare'](ps)
  expected=pd.read_csv(ROOT/'results/pdac/v3_phase5/program_comparisons.csv')
  got=rep[(rep.scope=='combined_metastasis')&(rep.program=='Hypoxia')].iloc[0];exp=expected[(expected.scope=='combined_metastasis')&(expected.program=='Hypoxia')].iloc[0]
  for key in ['n_pairs','p','q_bh','median_delta_metastasis_minus_primary']:
   reference.append({'field':key,'reproduced':float(got[key]),'saved_v3':float(exp[key]),'matches':bool(np.isclose(got[key],exp[key],atol=1e-12,rtol=0))})
  write_csv(pd.DataFrame(reference),out/'baseline_reproduction.csv',phase,'Recomputed v3 pipeline patient aggregation, exact test and five-program BH')
  if not all(r['matches'] for r in reference):raise RuntimeError('REPRODUCTION_MISMATCH saved pipeline patient test')
  errors=[]
  # Finish all raw-matrix reproduction checks before any sensitivity computation.
  for s in samples.itertuples(index=False):
   ad,counts=load(s);ad=normalize(ad[counts>=500].copy());values=score(ad,HYPO)
   saved=pd.read_csv(ROOT/'results/pdac/v3_phase5/slides'/f'{s.gsm}_spot_scores.csv.gz').set_index('barcode')
   target=saved.loc[ad.obs_names,'Hypoxia'].to_numpy();err=float(np.max(np.abs(values-target)))
   errors.append({'gsm':s.gsm,'n_spots':ad.n_obs,'max_abs_error':err,'matches':err<=1e-6})
   if err>1e-6:raise RuntimeError('REPRODUCTION_MISMATCH raw-matrix '+s.gsm)
   log.info('Raw reproduction %s max error %.3g',s.gsm,err);del ad;gc.collect()
  write_csv(pd.DataFrame(errors),out/'raw_reproduction.csv',phase,'Rescore every included raw slide; compare all baseline spot scores')
  html=(DATA/'sources/hallmark_hypoxia.html').read_text(encoding='iso-8859-1')
  hallmark=list(dict.fromkeys(re.findall(r"<a[^>]*href='[^']*ensembl[^']*'[^>]*>([^<]+)</a>",html)))
  if len(hallmark)!=200:raise RuntimeError('Independent Hallmark membership not verified as 200 displayed genes')
  write_json({'genes':hallmark,'source':'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_HYPOXIA.html','source_sha256':sha(DATA/'sources/hallmark_hypoxia.html'),'publication':'PMID 26771021','overlap_with_original':sorted(set(hallmark)&set(HYPO))},out/'independent_hypoxia_set.json',phase,'Opened public MSigDB member table and exact parsed membership')
  baseline=base[base.program=='Hypoxia'].copy();all_results=[];all_delta=[];slides=[];drops=[];means=[]
  r,d=paired(baseline,'baseline');all_results.append(r);all_delta.append(d)
  keep=samples.loc[samples.chemistry!='CytAssist_inferred','gsm'];r,d=paired(baseline[baseline.gsm.isin(keep)],'a_no_CytAssist');all_results.append(r);all_delta.append(d)
  variants=[('d_LOO_'+g,[x for x in HYPO if x!=g]) for g in HYPO]+[('d_Hallmark',hallmark)]
  for s in samples.itertuples(index=False):
   ad,counts=load(s);raw=ad.copy();ad=normalize(ad[counts>=500].copy())
   original=pd.read_csv(ROOT/'results/pdac/v3_phase5/slides'/f'{s.gsm}_spot_scores.csv.gz').set_index('barcode').loc[ad.obs_names,'Hypoxia'].to_numpy()
   common={'gsm':s.gsm,'patient':s.patient,'site':s.site}
   epithelial=[g for g in ['KRT19','KRT7','EPCAM'] if g in ad.var_names];host=[g for g in HOST[s.site] if g in ad.var_names]
   if len(epithelial)>=2 and len(host)>=2:
    e=np.asarray(ad[:,epithelial].X.mean(axis=1)).ravel();h=np.asarray(ad[:,host].X.mean(axis=1)).ravel()
    selected=(e>=np.quantile(e,2/3))&~((h>np.quantile(h,.9))&(h>0))
    if selected.sum()>=50:slides.append({**common,'scenario':'b_tumour_host','score':float(np.median(original[selected]))})
    else:drops.append({**common,'scenario':'b_tumour_host','reason':'fewer than 50 selected spots'})
    write_json({'gsm':s.gsm,'n_qc_spots':ad.n_obs,'n_selected':int(selected.sum()),'epithelial_genes':epithelial,'host_genes':host,'epithelial_threshold':float(np.quantile(e,2/3)),'host_threshold':float(np.quantile(h,.9))},out/'selection'/f'{s.gsm}.json',phase,'Tumour-rich and host-marker selection, declared before hypoxia examination')
   else:drops.append({**common,'scenario':'b_tumour_host','reason':'insufficient declared marker genes'})
   if np.mean(counts<2000)>.5:drops.append({**common,'scenario':'c_2000_UMI','reason':'more than half original spots below 2000 UMIs'})
   else:
    thin=raw[counts>=2000].copy();sc.pp.downsample_counts(thin,counts_per_cell=2000,random_state=42,replace=False)
    assert np.all(np.asarray(thin.X.sum(axis=1)).ravel()==2000)
    normalize(thin);slides.append({**common,'scenario':'c_2000_UMI','score':float(np.median(score(thin,HYPO)))})
    del thin
   for name,genes in variants:slides.append({**common,'scenario':name,'score':float(np.median(score(ad,genes)))})
   means.append(pd.Series(np.asarray(ad.X.sum(axis=0,dtype=np.float64)).ravel()/ad.n_obs,index=ad.var_names,name=s.gsm))
   log.info('Sensitivity scores %s complete',s.gsm);del ad,raw;gc.collect()
  tab=pd.DataFrame(slides)
  for scenario,g in tab.groupby('scenario',sort=False):r,d=paired(g,scenario);all_results.append(r);all_delta.append(d)
  write_csv(tab,out/'sensitivity_slide_scores.csv',phase,'Declared sensitivities; medians per slide')
  write_csv(pd.DataFrame(drops),out/'dropped_slides.csv',phase,'Every dropped sensitivity slide and reason')
  # Organ differences are descriptive and retain patient as the unit.
  site=baseline.groupby(['patient','site'],as_index=False).score.mean();org=[]
  for organ in sorted(site.site.unique()):
   if organ==PRIMARY:continue
   wide=site[site.site.isin([PRIMARY,organ])].pivot(index='patient',columns='site',values='score').dropna()
   for pat,row in wide.iterrows():org.append({'organ':organ,'patient':pat,'difference':row[organ]-row[PRIMARY]})
  org=pd.DataFrame(org);write_csv(org,out/'organ_paired_differences.csv',phase,'Organ-stratified descriptive paired differences; no inferential p-values')
  write_csv(org.groupby('organ',as_index=False).agg(n_pairs=('patient','size'),median_difference=('difference','median')),out/'organ_summary.csv',phase,'Descriptive organ sample sizes and effects')
  gene_means=pd.concat(means,axis=1,join='inner').mean(axis=1);write_csv(gene_means.rename('equal_slide_mean').rename_axis('gene').reset_index(),out/'gene_expression_matching.csv',phase,'Equal-slide expression used before random-set sampling')
  universe=gene_means[(gene_means>.01)&~gene_means.index.isin(HYPO)].sort_index()
  allbin=gene_means[gene_means>.01].sort_index();bins=np.minimum((allbin.rank(method='first')-1)*25/len(allbin),24).astype(int)
  rng=np.random.default_rng(42);randomsets=[];seen=set()
  while len(randomsets)<1000:
   chosen=[]
   for gene in HYPO:
    pool=[g for g in universe.index[bins.reindex(universe.index).to_numpy()==bins.loc[gene]] if g not in chosen]
    chosen.append(str(rng.choice(pool)))
   key=tuple(sorted(chosen))
   if key not in seen:seen.add(key);randomsets.append(chosen)
  write_json({'seed':42,'n_universe':len(universe),'sets':randomsets},out/'matched_random_gene_sets.json',phase,'1000 distinct size/mean-expression-matched signatures')
  randomslides=[]
  for i,s in enumerate(samples.itertuples(index=False)):
   ad,counts=load(s);ad=normalize(ad[counts>=500].copy());coef=coefficients(ad,randomsets)
   vals=np.asarray((ad.X@coef).toarray())
   # Check exact control-selection equivalence on actual data, not an invented fixture.
   if i==0:
    for k in [0,1,999]:
     expected=score(ad,randomsets[k],key='control_check')
     if not np.allclose(vals[:,k],expected,atol=1e-6,rtol=1e-6):raise ValueError('Random-set coefficient score disagrees with Scanpy')
   randomslides.append({'gsm':s.gsm,'patient':s.patient,'site':s.site,**{f'r{k}':v for k,v in enumerate(np.median(vals,axis=0))}})
   log.info('Random controls %s complete',s.gsm);del ad,coef,vals;gc.collect()
  rand=pd.DataFrame(randomslides);valuecols=[f'r{k}' for k in range(1000)]
  rsite=rand.groupby(['patient','site'])[valuecols].mean().reset_index();rsite['arm']=np.where(rsite.site==PRIMARY,'primary','metastasis')
  rarms=rsite.groupby(['patient','arm'])[valuecols].mean();rp=rarms.xs('primary',level='arm');rm=rarms.xs('metastasis',level='arm');effects=(rm-rp).median(axis=0).to_numpy()
  write_csv(pd.DataFrame({'random_set':range(1000),'median_paired_effect':effects}),out/'random_set_effects.csv',phase,'1000 patient-level paired random-signature effects')
  observed=all_results[0]['median_paired_effect'];tail=float((1+(effects>=observed).sum())/1001);percentile=float(100*np.mean(effects<=observed))
  tests=pd.DataFrame(all_results);tests['secondary_q_abc']=np.nan;abc=tests.scenario.isin(['a_no_CytAssist','b_tumour_host','c_2000_UMI']);tests.loc[abc,'secondary_q_abc']=bh(tests.loc[abc,'p'])
  write_csv(tests,out/'robustness_tests.csv',phase,'All prespecified paired tests and fixed a/b/c BH family')
  write_csv(pd.concat(all_delta,ignore_index=True),out/'patient_paired_differences.csv',phase,'Patient-level differences underlying every paired test')
  required=tests[tests.scenario.str.startswith(('a_','b_','c_','d_'))];failed=required[(required.median_paired_effect<=0)|(required.scenario.str.startswith(('a_','b_','c_'))&(required.p>=.05))]
  if (required.n_pairs<5).any():decision='INCONCLUSIVE'
  else:decision='ROBUST' if failed.empty else 'NOT ROBUST'
  fig,axes=plt.subplots(1,3,figsize=(15,4.5));plot=tests[~tests.scenario.str.startswith('d_LOO')]
  axes[0].barh(plot.scenario,plot.median_paired_effect,color='#326d8b');axes[0].axvline(0,color='black',lw=.6);axes[0].set_xlabel('Median paired score difference')
  loo=tests[tests.scenario.str.startswith('d_LOO')];axes[1].barh(loo.scenario.str.replace('d_LOO_',''),loo.median_paired_effect,color='#438d72');axes[1].axvline(0,color='black',lw=.6);axes[1].set_title('Leave-one-gene-out')
  axes[2].hist(effects,bins=40,color='#92adb9');axes[2].axvline(observed,color='#b94e42',label=f'Hypoxia; tail p={tail:.3g}');axes[2].legend();axes[2].set_title('1,000 matched random signatures');axes[2].set_xlabel('Median paired effect')
  fig.suptitle('Visium hypoxia robustness: '+decision);fig.tight_layout();figpaths=figure(fig,phase,'hypoxia_robustness');plt.close(fig)
  summary={'status':'DONE','decision':decision,'failed_checks':failed.scenario.tolist(),'baseline_reproduced':True,'n_slides':len(samples),'n_original_pairs':int(got.n_pairs),'random_tail_p':tail,'random_percentile':percentile,'n_random_sets':1000,'figure':figpaths,'results':tests.to_dict(orient='records'),'plan_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()}
  write_json(summary,out/'summary.json',phase,'Prespecified robustness verdict including negative/inconclusive checks')
  log_progress('## Phase 1 DONE\n\n'+decision+'. Failed checks: '+str(summary['failed_checks'])+'. Source: '+rel(out/'summary.json')+'. Baseline aggregate/test and every raw slide reproduced. Random signature tail p='+str(tail)+'.')
  log.info('Phase 1 %s; failed=%s; random tail=%s',decision,summary['failed_checks'],tail)
 except Exception as e:
  log.exception('Phase 1 failed')
  write_json({'status':'BLOCKED','reason':str(e)},out/'blocked.json',phase,'Phase 1 blocker; no invented robustness verdict')
  log_progress('## Phase 1 BLOCKED\n\n'+str(e)+'. Evidence: '+rel(out/'blocked.json'))
  raise
 finally:finish(log)
if __name__=='__main__':main()
