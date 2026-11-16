from common import *
import warnings
from scipy.stats import chi2,rankdata,spearmanr
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import SurvfuncRight

PH_SOURCE='https://raw.githubusercontent.com/CamDavidsonPilon/lifelines/master/lifelines/statistics.py'

def ph_test(fit,time,event,names):
 event=np.asarray(event)==1;res=np.asarray(fit.schoenfeld_residuals)[event];t=np.asarray(time)[event]
 if not np.isfinite(res).all():raise ValueError('Nonfinite event Schoenfeld residuals')
 n=len(t);scaled=n*(res@np.asarray(fit.cov_params()));rank=rankdata(t,method='average');rank=rank-rank.mean()
 stat=(np.sum(rank[:,None]*scaled,axis=0)**2)/(n*np.asarray(fit.bse)**2*np.sum(rank**2));p=chi2.sf(stat,1)
 table=pd.DataFrame({'term':names,'statistic':stat,'p':p,'n_events':n,'diagnostic':'rank-time scaled Schoenfeld score approximation; average tied ranks','source':PH_SOURCE})
 residual=pd.DataFrame(scaled,columns=['residual_'+name for name in names]);residual['time_months']=t;residual['event_rank']=rank
 return table,residual

def validate_ph():
 rng=np.random.default_rng(42);records=[]
 for mode in ['proportional','time_varying']:
  for rep in range(100):
   n=300;z=rng.binomial(1,.5,n);u=rng.exponential(1,n);early=.04*np.exp(.6*z if mode=='proportional' else 1.4*z);late=.04*np.exp(.6*z if mode=='proportional' else -1.4*z);t0=8
   latent=np.where(u<=early*t0,u/early,t0+(u-early*t0)/late);c=rng.exponential(30,n);t=np.minimum(latent,c);e=(latent<=c).astype(int)
   fit=PHReg(t,z[:,None],status=e,ties='breslow').fit(disp=False,maxiter=100)
   table,_=ph_test(fit,t,e,['z'])
   step=1e-5;b=fit.params;finite=(fit.model.loglike(b+step)-fit.model.loglike(b-step))/(2*step);grad=fit.model.score(b)[0]
   records.append({'mode':mode,'replicate':rep,'p':float(table.p.iloc[0]),'beta':float(b[0]),'gradient_error':abs(finite-grad)})
 result=pd.DataFrame(records);null=float((result[result['mode']=='proportional'].p<.05).mean());power=float((result[result['mode']=='time_varying'].p<.05).mean())
 passed=null<=.15 and power>=.5 and result.gradient_error.max()<1e-4
 return result,{'status':'PASS' if passed else 'BLOCKED','null_rejection_rate':null,'time_varying_detection_rate':power,'max_gradient_error':float(result.gradient_error.max()),'source':PH_SOURCE,'limitations':'Approximate score diagnostic; not exact equivalence to current R cox.zph, especially for ties'}

def score_signatures(expression,gene_sets):
 # Patient rows, gene columns; callers preserve source scale and cohort identity.
 out=pd.DataFrame(index=expression.index);coverage=[];components={}
 for name,genes in gene_sets.items():
  present=[g for g in genes if g in expression.columns and np.isfinite(expression[g]).all() and expression[g].std(ddof=1)>0]
  valid=len(present)>=5 and len(present)>=.5*len(genes)
  coverage.append({'signature':name,'n_declared':len(genes),'n_present_finite_variable':len(present),'coverage_fraction':len(present)/len(genes),'genes_used':';'.join(present),'missing_or_nonfinite_or_constant':';'.join(g for g in genes if g not in present),'status':'SCORABLE' if valid else 'BLOCKED_coverage'})
  if valid:
   values=expression[present];z=(values-values.mean())/values.std(ddof=1);raw=z.mean(axis=1);components[name]=raw
   out[name]=(raw-raw.mean())/raw.std(ddof=1)
 if 'CAF' in components and 'T_cell_proxy' in components:
  raw=components['CAF']-components['T_cell_proxy'];out['Immune_exclusion_proxy']=(raw-raw.mean())/raw.std(ddof=1)
  coverage.append({'signature':'Immune_exclusion_proxy','n_declared':len(gene_sets['CAF'])+len(gene_sets['T_cell_proxy']),'n_present_finite_variable':sum(r['n_present_finite_variable'] for r in coverage if r['signature'] in ['CAF','T_cell_proxy']),'coverage_fraction':sum(r['n_present_finite_variable'] for r in coverage if r['signature'] in ['CAF','T_cell_proxy'])/(len(gene_sets['CAF'])+len(gene_sets['T_cell_proxy'])),'genes_used':'CAF mean gene z minus T-cell mean gene z','missing_or_nonfinite_or_constant':'See component coverage','status':'SCORABLE'})
 return out,pd.DataFrame(coverage)

def fit_models(clinical,scores,signatures,cohort):
 result=[];coef=[];ph=[];residuals=[];membership=[]
 d=clinical.set_index('PATIENT_ID').join(scores,validate='one_to_one')
 for signature in signatures:
  for family in ['unadjusted','clinical','clinical_purity']:
   required=[signature,'time_months','event']+(['age','sex','stage'] if family!='unadjusted' else [])+(['purity'] if family=='clinical_purity' else [])
   row={'cohort':cohort,'signature':signature,'model':family,'n_planned_signatures':len(signatures),'status':'BLOCKED','reason':''}
   if signature not in d:
    row['reason']='insufficient signature gene coverage';result.append(row);continue
   selected=d.dropna(subset=required).copy();row['n']=len(selected);row['events']=int(selected.event.sum());row['excluded_complete_case']=len(d)-len(selected)
   x=pd.DataFrame({'signature':selected[signature]},index=selected.index)
   if family!='unadjusted':
    x['age_decades']=selected.age/10;x['male']=selected.sex.eq('Male').astype(int);x['stage_II']=selected.stage.eq('II').astype(int);x['stage_III_IV']=selected.stage.eq('III_IV').astype(int)
   if family=='clinical_purity':x['purity']=selected.purity
   row['n_parameters']=x.shape[1];row['events_per_parameter']=row['events']/x.shape[1]
   membership.extend({'cohort':cohort,'signature':signature,'model':family,'patient':patient,'included':patient in selected.index} for patient in d.index)
   if row['events']<10*x.shape[1]:row['reason']='fewer than 10 events per fitted coefficient';result.append(row);continue
   if np.linalg.matrix_rank(x.to_numpy())!=x.shape[1]:row['reason']='rank-deficient fixed model';result.append(row);continue
   try:
    with warnings.catch_warnings(record=True) as warning_list:
     warnings.simplefilter('always');fit=PHReg(selected.time_months.to_numpy(),x.to_numpy(),status=selected.event.to_numpy(),ties='breslow').fit(disp=False,maxiter=100)
    row['warnings']='; '.join(str(w.message) for w in warning_list)
    if not np.isfinite(fit.params).all() or not np.isfinite(fit.bse).all() or np.linalg.norm(fit.model.score(fit.params))>1e-4:raise ValueError('Nonfinite or nonconverged Cox fit')
    row.update({'status':'DONE','hr':float(np.exp(fit.params[0])),'ci_low':float(np.exp(fit.params[0]-1.95996398454*fit.bse[0])),'ci_high':float(np.exp(fit.params[0]+1.95996398454*fit.bse[0])),'p':float(fit.pvalues[0]),'beta':float(fit.params[0]),'gradient_norm':float(np.linalg.norm(fit.model.score(fit.params)))})
    for k,term in enumerate(x.columns):coef.append({'cohort':cohort,'signature':signature,'model':family,'term':term,'beta':float(fit.params[k]),'se':float(fit.bse[k]),'p':float(fit.pvalues[k])})
    pt,res=ph_test(fit,selected.time_months.to_numpy(),selected.event.to_numpy(),list(x.columns))
    for key,value in [('cohort',cohort),('signature',signature),('model',family)]:pt[key]=value;res[key]=value
    ph.append(pt);residuals.append(res)
   except Exception as e:row['status']='BLOCKED';row['reason']=str(e)
   result.append(row)
 tab=pd.DataFrame(result)
 if 'p' not in tab:tab['p']=np.nan
 tab['q_bh']=np.nan;tab['n_tested_in_family']=0
 for _,ix in tab.groupby(['cohort','model']).groups.items():
  pv=tab.loc[ix,'p'];q=bh(pv.fillna(1).to_numpy());q[~pv.notna().to_numpy()]=np.nan
  tab.loc[ix,'q_bh']=q;tab.loc[ix,'n_tested_in_family']=int(pv.notna().sum())
 pt=pd.concat(ph,ignore_index=True) if ph else pd.DataFrame()
 if len(pt):pt['q_bh_PH']=pt.groupby(['cohort','model']).p.transform(lambda v:pd.Series(bh(v.to_numpy()),index=v.index))
 return tab,pd.DataFrame(coef),pt,pd.concat(residuals,ignore_index=True) if residuals else pd.DataFrame(),pd.DataFrame(membership)

def km_data(clinical,scores,signatures,cohort):
 d=clinical.set_index('PATIENT_ID').join(scores,validate='one_to_one');steps=[];risk=[];cuts=[]
 for name in signatures:
  if name not in d:continue
  cut=float(d[name].median());cuts.append({'cohort':cohort,'signature':name,'median_cutoff':cut})
  for label,mask in [('low',d[name]<=cut),('high',d[name]>cut)]:
   g=d[mask];fit=SurvfuncRight(g.time_months.to_numpy(),g.event.to_numpy())
   steps.append({'cohort':cohort,'signature':name,'group':label,'time':0.,'survival':1.,'n_group':len(g)})
   for t,p in zip(fit.surv_times,fit.surv_prob):steps.append({'cohort':cohort,'signature':name,'group':label,'time':float(t),'survival':float(p),'n_group':len(g)})
   for t in [0,12,24,36,60]:risk.append({'cohort':cohort,'signature':name,'group':label,'time':t,'at_risk':int((g.time_months>=t).sum())})
 return pd.DataFrame(steps),pd.DataFrame(risk),pd.DataFrame(cuts)

def correlations(clinical,scores,signatures,cohort):
 d=clinical.set_index('PATIENT_ID').join(scores);cols=[s for s in signatures if s in scores]+['purity'];records=[]
 for i,a in enumerate(cols):
  for b in cols[i+1:]:
   g=d[[a,b]].dropna();rho,p=spearmanr(g[a],g[b]) if len(g)>2 else (np.nan,np.nan)
   records.append({'cohort':cohort,'signature_a':a,'signature_b':b,'n':len(g),'spearman_rho':float(rho),'p_descriptive':float(p),'mostly_purity_proxy_by_declared_rule':b=='purity' and abs(rho)>=.7})
 return pd.DataFrame(records)
