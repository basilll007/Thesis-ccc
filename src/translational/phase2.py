from common import *
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.duration.survfunc import SurvfuncRight

def source_file(cohort,name):
 p=DATA/cohort/'lfs'/name
 return p if p.exists() else DATA/cohort/name

def stage_group(s):
 if not isinstance(s,str):return np.nan
 x=s.upper().replace('STAGE','').strip()
 if x.startswith(('III','IV')):return 'III_IV'
 if x.startswith('II'):return 'II'
 if x.startswith('I'):return 'I'
 return np.nan

def cellmean(s):
 try:
  values=[float(x) for x in str(s).split(';')]
  if not values or not all(0<=x<=100 for x in values):return np.nan
  return float(np.mean(values))/100
 except (ValueError,TypeError):return np.nan

def main():
 out,log=start('phase2');phase='phase2'
 try:
  sources=[]
  for name in ['extra_downloads.json','cptac_downloads.json']:
   if (DATA/name).exists():sources.extend(json.loads((DATA/name).read_text()))
  acq=sorted((ROOT/'results/translational/acquisition').glob('*/downloads.json'))[-1];sources.extend(json.loads(acq.read_text()))
  write_csv(pd.DataFrame(sources),out/'source_downloads.csv',phase,'All actual source/download identities and CAPTCHA block; manual URL instructions retained')
  pur=pd.read_csv(DATA/'gdc/TCGA_mastercalls.abs_tables_JSedit.fixed.txt',sep='\t')
  pur=pur[(pur['call status']=='called')&pur.purity.between(0,1)].groupby('array').purity.mean()
  cohorts=[];summaries=[];missing=[];vet=[];config=[];exclusions=[]
  for name,dirn in [('TCGA','paad_pancan'),('CPTAC','paad_cptac')]:
   p=source_file(dirn,'data_clinical_patient.txt');s=source_file(dirn,'data_clinical_sample.txt');e=source_file(dirn,'data_mrna_seq_v2_rsem.txt')
   pat=pd.read_csv(p,sep='\t',comment='#');sam=pd.read_csv(s,sep='\t',comment='#')
   assert pat.PATIENT_ID.is_unique and sam.SAMPLE_ID.is_unique
   exprcols=pd.read_csv(e,sep='\t',nrows=0).columns.tolist();expression_samples=[c for c in exprcols if c not in ['Hugo_Symbol','Entrez_Gene_Id']]
   merged=sam.merge(pat,on='PATIENT_ID',validate='many_to_one',suffixes=('','_patient'))
   merged['has_expression']=merged.SAMPLE_ID.isin(expression_samples)
   merged=merged.sort_values('SAMPLE_ID');merged['one_sample']=~merged.PATIENT_ID.duplicated()
   if name=='TCGA':
    merged['strict_pdac']=merged.ICD_O_3_HISTOLOGY.eq('8500/3');merged['histology_evidence']=merged.ICD_O_3_HISTOLOGY
    merged['time_months']=pd.to_numeric(merged.OS_MONTHS,errors='coerce');merged['event']=merged.OS_STATUS.map({'1:DECEASED':1,'0:LIVING':0})
    merged['purity']=merged.SAMPLE_ID.map(pur);merged['purity_method']='GDC ABSOLUTE called aliquots, mean per array'
    merged['stage']=merged.AJCC_PATHOLOGIC_TUMOR_STAGE.map(stage_group);merged['grade']=merged.GRADE
    merged['is_primary']=merged.SAMPLE_TYPE.eq('Primary');exprscale='linear_RSEM_log2p1'
    citation='TCGA PanCancer Atlas Cell 2018; TCGA Research Network 2017 doi:10.1016/j.ccell.2017.07.007; Cerami 2012 and Gao 2013 cBioPortal'
    hist_source='ICD_O_3_HISTOLOGY=8500/3; NCI morphology table; published curated sample list unavailable (CAPTCHA)'
   else:
    merged['strict_pdac']=merged.HISTOLOGY_DIAGNOSIS.eq('Pancreatic Ductal Adenocarcinoma');merged['histology_evidence']=merged.HISTOLOGY_DIAGNOSIS
    merged['time_months']=pd.to_numeric(merged.FOLLOW_UP_DAYS,errors='coerce')/(365.25/12);merged['event']=merged.VITAL_STATUS.map({'Deceased':1,'Living':0})
    merged['purity']=merged.NEOPLASTIC_CELLULARITY.map(cellmean);merged['purity_method']='Pathology cellularity mean of semicolon values /100; percentage scale inferred'
    merged['stage']=merged.TUMOR_STAGE_PATHOLOGICAL.map(stage_group);merged['grade']=np.nan
    merged['is_primary']=True;exprscale='already_log2_RSEM_UQ'
    citation='Cao et al. Cell 2021 doi:10.1016/j.cell.2021.08.023 PMID34534465; Cerami 2012 and Gao 2013 cBioPortal'
    hist_source='HISTOLOGY_DIAGNOSIS explicitly Pancreatic Ductal Adenocarcinoma; excludes adenosquamous'
   merged['age']=pd.to_numeric(merged.AGE,errors='coerce');merged['sex']=merged.SEX.map({'Male':'Male','Female':'Female'})
   merged['eligible']=merged.has_expression&merged.one_sample&merged.is_primary&merged.time_months.gt(0)&merged.event.isin([0,1])
   for row in merged.itertuples(index=False):
    if not row.eligible:exclusions.append({'dataset':name,'patient':row.PATIENT_ID,'sample':row.SAMPLE_ID,'expression':row.has_expression,'one_sample':row.one_sample,'primary':row.is_primary,'time_months':row.time_months,'event':row.event})
   merged['dataset']=name
   keepcols=['dataset','PATIENT_ID','SAMPLE_ID','strict_pdac','histology_evidence','has_expression','is_primary','time_months','event','age','sex','stage','grade','purity','purity_method','eligible']
   write_csv(merged[keepcols],out/f'{name}_cohort_all_records.csv',phase,'Source clinical records, explicit histology and every eligibility decision')
   for cohort,mask in [('all_primary',merged.eligible),('confirmed_PDAC',merged.eligible&merged.strict_pdac)]:
    d=merged[mask].copy();cohort_id=name+'_'+cohort
    assert d.PATIENT_ID.is_unique and d.SAMPLE_ID.is_unique
    follow=SurvfuncRight(d.time_months.to_numpy(),1-d.event.to_numpy())
    at=follow.surv_times[follow.surv_prob<=.5];median_follow=float(at[0]) if len(at) else np.nan
    summaries.append({'cohort':cohort_id,'n_patients':len(d),'n_events':int(d.event.sum()),'n_censored':int((d.event==0).sum()),'median_observed_time_months':float(d.time_months.median()),'reverse_KM_median_followup_months':median_follow,'n_with_purity':int(d.purity.notna().sum()),'n_source_patients':len(pat),'n_expression_samples':len(expression_samples)})
    for col in ['time_months','event','age','sex','stage','grade','purity']:missing.append({'cohort':cohort_id,'field':col,'n_missing':int(d[col].isna().sum()),'fraction_missing':float(d[col].isna().mean())})
    q=out/f'{cohort_id}_patients.csv';write_csv(d[keepcols],q,phase,'One expression-profiled patient per row; eligibility and histology declared before survival models')
    config.append({'cohort':cohort_id,'dataset':name,'clinical_table':rel(q),'expression_file':str(e),'expression_sha256':sha(e),'expression_scale':exprscale,'strict':cohort=='confirmed_PDAC'})
   license_text=(DATA/dirn/'LICENSE').read_text(encoding='utf-8')
   rawbase='https://raw.githubusercontent.com/cBioPortal/datahub/master/public/'+('paad_tcga_pan_can_atlas_2018' if name=='TCGA' else 'paad_cptac_2021')+'/'
   vet.append({'dataset':name,'status':'OBTAINED','n_source_patients':len(pat),'n_source_samples':len(sam),'n_expression_samples':len(expression_samples),'data_type':exprscale,'clinical_fields':';'.join(pat.columns)+';'+(';'.join(sam.columns)),'histology_restriction':hist_source,'access_terms':'Study licence references Broad GDAC and preliminary-data caution' if name=='TCGA' else 'ODbL: attribution, open derivatives and share-alike','terms_url':rawbase+'LICENSE','metadata_url':rawbase+'meta_study.txt','expression_metadata_url':rawbase+'meta_mrna_seq_v2_rsem.txt','citation_requirements':citation,'purity_source':'https://gdc.cancer.gov/about-data/publications/pancanatlas' if name=='TCGA' else rawbase+'data_clinical_sample.txt','licence_sha256':sha(DATA/dirn/'LICENSE')})
   log.info('%s: source patients=%d, expression samples=%d',name,len(pat),len(expression_samples))
  write_csv(pd.DataFrame(vet),ROOT/'results/translational/cohort_vetting.csv',phase,'Verified source counts, fields, access terms, histology restrictions and citation requirements')
  write_csv(pd.DataFrame(summaries),out/'cohort_summary.csv',phase,'Patient/event counts, reverse-KM follow-up and observed-time medians')
  write_csv(pd.DataFrame(missing),out/'missingness.csv',phase,'Covariate missingness per exact cohort')
  write_csv(pd.DataFrame(exclusions),out/'exclusions.csv',phase,'Excluded clinical records; expression and endpoint mismatches remain explicit')
  write_json(config,out/'cohort_config.json',phase,'Immutable clinical and expression identities for Phase 3')
  fig,axes=plt.subplots(1,2,figsize=(13,4.5));ss=pd.DataFrame(summaries)
  axes[0].barh(ss.cohort,ss.n_patients,color='#74a0b1',label='Patients');axes[0].barh(ss.cohort,ss.n_events,color='#bd6652',label='Events');axes[0].legend();axes[0].set_title('Verified analysis cohorts')
  mm=pd.DataFrame(missing).pivot(index='cohort',columns='field',values='fraction_missing');im=axes[1].imshow(mm,aspect='auto',vmin=0,vmax=1,cmap='Blues');axes[1].set_xticks(range(len(mm.columns)),mm.columns,rotation=45);axes[1].set_yticks(range(len(mm)),mm.index);axes[1].set_title('Missing covariates');fig.colorbar(im,ax=axes[1],label='Fraction missing');fig.tight_layout();figpaths=figure(fig,phase,'cohort_vetting');plt.close(fig)
  write_json({'status':'DONE','n_datasets':len(vet),'cohorts':summaries,'config':rel(out/'cohort_config.json'),'figure':figpaths,'curated_TCGA_supplement':'BLOCKED (CAPTCHA); clinical histology fallback used','CPTAC_limitations':['Pathology cellularity differs from ABSOLUTE; percent scale inferred','FOLLOW_UP_DAYS and vital status used for OS; source paper 30-day exclusion not applied','Treatment and grade not modelled']},out/'summary.json',phase,'Phase 2 completion with source and eligibility limitations')
  log_progress('## Phase 2 DONE\n\nTwo public cohorts acquired. Published TCGA curated supplement is BLOCKED by CAPTCHA; strict clinical histology fallback and all-primary sensitivity used. Cohort counts and follow-up: '+rel(out/'cohort_summary.csv')+'.')
 finally:finish(log)
if __name__=='__main__':main()
