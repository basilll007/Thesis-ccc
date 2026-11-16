from common import *
import ast,gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from survival import score_signatures,fit_models,km_data,correlations,validate_ph

SIGNATURES=['Hypoxia','CAF','EMT','Interferon','Immune_infiltration','Immune_exclusion_proxy']

def gene_sets():
 tree=ast.parse((ROOT/'src/pdac/v3_phase5_programs.py').read_text(encoding='utf-8'));ns={}
 for node in tree.body:
  if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='PROGRAMS' for t in node.targets):exec(compile(ast.Module(body=[node],type_ignores=[]),'v3_programs','exec'),ns)
 sets=ns['PROGRAMS'];sets['T_cell_proxy']='CD3D CD3E CD3G CD2 TRAC CD8A CD4 IL7R CD247'.split();return sets

def load_expression(conf,clinical):
 p=Path(conf['expression_file'])
 if sha(p)!=conf['expression_sha256']:raise ValueError('Expression input identity changed')
 columns=['Hugo_Symbol']+clinical.SAMPLE_ID.tolist();raw=pd.read_csv(p,sep='\t',usecols=columns)
 raw=raw.dropna(subset=['Hugo_Symbol']);raw['Hugo_Symbol']=raw.Hugo_Symbol.astype(str)
 values=raw.groupby('Hugo_Symbol',sort=True)[clinical.SAMPLE_ID.tolist()].mean().T
 if conf['expression_scale']=='linear_RSEM_log2p1':
  if np.nanmin(values.to_numpy())<0:raise ValueError('Negative values in declared linear RSEM source')
  values=np.log2(values+1)
 elif conf['expression_scale']!='already_log2_RSEM_UQ':raise ValueError('Unknown expression scale')
 values.index=clinical.set_index('SAMPLE_ID').loc[values.index,'PATIENT_ID'].to_numpy()
 assert values.index.is_unique and set(values.index)==set(clinical.PATIENT_ID)
 return values

def cohort_configs():
 p=sorted((ROOT/'results/translational/phase2').glob('*/cohort_config.json'))[-1]
 return json.loads(p.read_text())

def main():
 out,log=start('phase3');phase='phase3'
 try:
  validation,receipt=validate_ph();write_csv(validation,out/'PH_synthetic_validation.csv',phase,'200 prespecified proportional/time-varying examples and finite-difference gradient validation')
  write_json(receipt,out/'PH_validation_receipt.json',phase,'Independent rank-time PH diagnostic validation and explicit approximation limits')
  if receipt['status']!='PASS':raise RuntimeError('PH diagnostic synthetic validation failed; no patient outcome models executed')
  sets=gene_sets();write_json(sets,out/'declared_gene_sets.json',phase,'Verbatim v3 programs and fixed T-cell proxy; no outcome-selected genes')
  models=[];coefficients=[];ph=[];resids=[];membership=[];corr=[];km=[];risk=[];cutoffs=[];coverage=[];scored={};config=cohort_configs()
  for conf in config:
   name=conf['cohort'];clinical=pd.read_csv(ROOT/conf['clinical_table']);expression=load_expression(conf,clinical)
   scores,cov=score_signatures(expression,sets);cov['cohort']=name;coverage.append(cov)
   write_csv(scores.rename_axis('PATIENT_ID').reset_index(),out/f'{name}_patient_scores.csv',phase,'Per-patient gene-z program scores, scaled per cohort; no optimised cutpoints')
   m,c,p,r,inc=fit_models(clinical,scores,SIGNATURES,name)
   models.append(m);coefficients.append(c);ph.append(p);resids.append(r);membership.append(inc)
   correlation=correlations(clinical,scores,SIGNATURES,name);corr.append(correlation)
   curves,at_risk,cuts=km_data(clinical,scores,SIGNATURES,name);km.append(curves);risk.append(at_risk);cutoffs.append(cuts)
   if conf['strict']:scored[conf['dataset']]={'cohort':name,'models':m,'correlation':correlation,'curves':curves,'clinical':clinical}
   log.info('%s: n=%d, model status %s',name,len(clinical),m.status.value_counts().to_dict());del expression;gc.collect()
  def cat(frames):return pd.concat([f for f in frames if len(f)],ignore_index=True) if any(len(f) for f in frames) else pd.DataFrame()
  m=cat(models);p=cat(ph);c=cat(corr)
  for df,name,desc in [(m,'cox_models.csv','Six signatures by three fixed Cox models; BH includes six planned tests per cohort/model'),(cat(coefficients),'cox_all_coefficients.csv','All fitted signature and clinical coefficients'),(p,'PH_diagnostics.csv','Approximate rank-time Schoenfeld diagnostics; BH across all PH tests per family'),(cat(resids),'scaled_schoenfeld_residuals.csv','Event-level scaled Schoenfeld residuals and transformed time'),(cat(membership),'model_patient_membership.csv','Complete-case membership for every patient/model'),(c,'signature_correlations.csv','Purity and inter-signature correlations and declared proxy flag'),(cat(coverage),'gene_coverage.csv','Declared, present, finite and variable genes per signature/cohort'),(cat(km),'KM_curves.csv','Descriptive median-split survival step curves'),(cat(risk),'KM_at_risk.csv','At-risk counts at fixed month intervals'),(cat(cutoffs),'KM_cutoffs.csv','Fixed per-cohort median cutpoints')]:write_csv(df,out/name,phase,desc)
  directions=[]
  tcga=m[m.cohort.eq('TCGA_confirmed_PDAC')];cptac=m[m.cohort.eq('CPTAC_confirmed_PDAC')]
  for signature in SIGNATURES:
   for family in ['unadjusted','clinical','clinical_purity']:
    a=tcga[(tcga.signature==signature)&(tcga.model==family)].iloc[0];b=cptac[(cptac.signature==signature)&(cptac.model==family)].iloc[0]
    valid=a.status=='DONE' and b.status=='DONE'
    directions.append({'signature':signature,'model':family,'TCGA_status':a.status,'CPTAC_status':b.status,'TCGA_HR':a.get('hr',np.nan),'CPTAC_HR':b.get('hr',np.nan),'same_direction':'yes' if valid and np.sign(a.beta)==np.sign(b.beta) else ('no' if valid else 'UNTESTABLE'),'TCGA_q':a.get('q_bh',np.nan),'CPTAC_q':b.get('q_bh',np.nan)})
  write_csv(pd.DataFrame(directions),out/'external_validation_directions.csv',phase,'Unchanged-model sign comparison in strict TCGA and CPTAC; no post-hoc tuning')
  fig=plt.figure(figsize=(21,11));grid=fig.add_gridspec(3,6,height_ratios=[1,1,1.25])
  for row,dataset in enumerate(['TCGA','CPTAC']):
   info=scored[dataset]
   for col,signature in enumerate(SIGNATURES):
    ax=fig.add_subplot(grid[row,col]);g=info['curves'][info['curves'].signature==signature]
    for label,color in [('low','#287697'),('high','#bd5946')]:
     v=g[g.group==label];ax.step(v.time,v.survival,where='post',color=color,lw=1.2,label=label+f' n={int(v.n_group.iloc[0])}')
    ax.set(title=dataset+' '+signature.replace('_',' '),xlabel='Months',ylim=(0,1.05));ax.legend(fontsize=6);ax.spines[['top','right']].set_visible(False)
    if col==0:ax.set_ylabel('Descriptive OS probability')
  ax=fig.add_subplot(grid[2,:3]);labels=[];counter=0
  for signature in SIGNATURES:
   for dataset,color,offset in [('TCGA','#287697',0),('CPTAC','#bd5946',.28)]:
    sub=m[(m.cohort==dataset+'_confirmed_PDAC')&(m.signature==signature)&(m.model=='clinical_purity')]
    y=counter+offset
    if len(sub) and sub.status.iloc[0]=='DONE':
     rr=sub.iloc[0];ax.errorbar(rr.hr,y,xerr=[[rr.hr-rr.ci_low],[rr.ci_high-rr.hr]],fmt='o',color=color,ms=4)
    labels.append((y,dataset+' '+signature));
   counter+=1
  ax.set_yticks([x[0] for x in labels],[x[1] for x in labels],fontsize=7);ax.axvline(1,color='black',ls='--',lw=.6);ax.set_xscale('log');ax.set_xlabel('HR per SD, age/sex/stage/purity adjusted');ax.set_title('Strict ductal cohorts; fixed models')
  ax=fig.add_subplot(grid[2,3:5]);co=scored['TCGA']['correlation'];names=SIGNATURES+['purity'];mat=pd.DataFrame(np.eye(len(names)),index=names,columns=names)
  for rr in co.itertuples():mat.loc[rr.signature_a,rr.signature_b]=mat.loc[rr.signature_b,rr.signature_a]=rr.spearman_rho
  im=ax.imshow(mat,vmin=-1,vmax=1,cmap='RdBu_r');ax.set_xticks(range(len(names)),names,rotation=90,fontsize=6);ax.set_yticks(range(len(names)),names,fontsize=6);ax.set_title('TCGA expression/purity correlations');fig.colorbar(im,ax=ax,fraction=.04)
  ax=fig.add_subplot(grid[2,5]);res=cat(resids);g=res[(res.cohort=='TCGA_confirmed_PDAC')&(res.signature=='Hypoxia')&(res.model=='clinical_purity')]
  if len(g):
   ax.scatter(g.time_months,g.residual_signature,s=5,alpha=.35)
  ax.set_title('Hypoxia PH diagnostic');ax.set_xlabel('Event month');ax.set_ylabel('Scaled Schoenfeld residual')
  fig.suptitle('Spatial program signatures in public bulk PDAC cohorts: expression associations',y=1.01);fig.tight_layout();figpaths=figure(fig,phase,'signature_survival');plt.close(fig)
  strict=m[m.cohort.eq('TCGA_confirmed_PDAC')&m.model.eq('clinical_purity')]
  significant=strict[(strict.status=='DONE')&(strict.q_bh<.05)]
  summary={'status':'DONE','n_primary_model_families':3,'n_signatures_per_family':6,'n_planned_primary_models':18,'n_all_cohort_models':len(m),'n_models_done':int((m.status=='DONE').sum()),'n_models_blocked':int((m.status!='DONE').sum()),'primary_adjusted_significant':significant.signature.tolist(),'primary_adjusted':strict.to_dict(orient='records'),'figure':figpaths,'PH_validation':receipt,'limitations':['Approximate PH diagnostic; report failures and do not establish PH from a nonsignificant test','CPTAC pathology cellularity differs from DNA ABSOLUTE purity','No treatment modelling or optimised cutpoints','Bulk expression cannot verify a spatial mechanism']}
  write_json(summary,out/'summary.json',phase,'Phase 3 fixed-model outcome, multiplicity and limitations')
  log_progress('## Phase 3 DONE\n\nFixed six-signature models across strict and all-primary TCGA/CPTAC cohorts. Primary purity-adjusted BH-significant signatures: '+str(summary['primary_adjusted_significant'])+'. Source: '+rel(out/'cox_models.csv')+'.')
 except Exception as e:
  log.exception('Phase 3 failure');write_json({'status':'BLOCKED','reason':str(e)},out/'blocked.json',phase,'Phase 3 execution or validity blocker');log_progress('## Phase 3 BLOCKED\n\n'+str(e)+'. '+rel(out/'blocked.json'));raise
 finally:finish(log)
if __name__=='__main__':main()
