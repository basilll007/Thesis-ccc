from common import *
import re

def main():
 out,log=start('acquisition');records=[]
 raw='https://raw.githubusercontent.com/cBioPortal/datahub/master/public/paad_tcga_pan_can_atlas_2018/'
 files=['meta_study.txt','LICENSE','README.md','meta_mrna_seq_v2_rsem.txt','data_clinical_patient.txt','data_clinical_sample.txt','data_mrna_seq_v2_rsem.txt']
 for name in files:
  p,rec=download('paad_pancan/'+name,raw+name,log);records.append(rec)
  if p and p.stat().st_size<300:
   content=p.read_text()
   if content.startswith('version https://git-lfs.github.com/spec/v1'):
    expected=re.search(r'oid sha256:(\w+)',content).group(1)
    target='https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/paad_tcga_pan_can_atlas_2018/'+name
    q,rr=download('paad_pancan/lfs/'+name,target,log);records.append(rr)
    if q and sha(q)!=expected:raise ValueError('Git LFS object SHA mismatch: '+name)
 for name,url in [
 ('sources/hallmark_hypoxia.html','https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_HYPOXIA.html'),
 ('sources/cbio_faq.html','https://docs.cbioportal.org/user-guide/faq/'),
 ('sources/tcga2017_paper.html','https://pmc.ncbi.nlm.nih.gov/articles/PMC5964983/')]:
  p,rec=download(name,url,log);records.append(rec)
 write_json(records,out/'downloads.json','acquisition','Every download URL, content identity or blocker; Git LFS identities verified')
 write_csv(pd.DataFrame(records),out/'downloads.csv','acquisition','Opened download/source URLs and blocks')
 log_progress('## Acquisition\n\nSource downloads and access blocks: '+rel(out/'downloads.csv')+'. Native Git LFS media endpoints are used only when a successfully read raw pointer specifies the exact expected SHA-256; no denied endpoint is bypassed.')
 finish(log)
if __name__=='__main__':main()
