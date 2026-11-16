from __future__ import annotations
import os,sys,csv,json,hashlib,logging,subprocess,urllib.request,urllib.error,importlib.metadata
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
SOURCE=Path('F:/Thesis')
DATA=SOURCE/'data/tcga'
CACHE=ROOT/'results/translational/cache'
CACHE.mkdir(parents=True,exist_ok=True)
sys.dont_write_bytecode=True
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='4'
os.environ['NUMBA_CACHE_DIR']=str(CACHE/'numba')
os.environ['MPLCONFIGDIR']=str(CACHE/'matplotlib')
if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
import numpy as np,pandas as pd

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def rel(p):return Path(p).resolve().relative_to(ROOT).as_posix()
def register(p,phase,description):
 p=Path(p);man=ROOT/'results/translational/MANIFEST.csv'
 old=[]
 if man.exists():
  with man.open(encoding='utf-8',newline='') as f:old=list(csv.DictReader(f))
 row={'path':rel(p),'phase':phase,'description':description,'bytes':p.stat().st_size,'sha256':sha(p),'registered_at':datetime.now(timezone.utc).isoformat()}
 old=[r for r in old if r['path']!=row['path']]+[row]
 with man.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(row));w.writeheader();w.writerows(old)
def start(phase):
 stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
 out=ROOT/f'results/translational/{phase}'/stamp;out.mkdir(parents=True)
 logger=logging.getLogger(phase+stamp);logger.setLevel(logging.INFO)
 file=out/'run.log';h=logging.FileHandler(file,encoding='utf-8');s=logging.StreamHandler(sys.stdout)
 for x in [h,s]:x.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'));logger.addHandler(x)
 logger.file=file;logger.phase=phase
 logger.info('commit=%s seed=42 max_threads=4',subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
 vers={}
 for name in ['numpy','pandas','scipy','scanpy','anndata','statsmodels','h5py','matplotlib','lifelines']:
  try:vers[name]=importlib.metadata.version(name)
  except importlib.metadata.PackageNotFoundError:vers[name]='not installed'
 logger.info('python=%s versions=%s',sys.version,vers)
 return out,logger
def finish(log):
 for h in log.handlers:h.flush()
 register(log.file,log.phase,'Run log: parameters, package versions, commit and lifecycle')
def log_progress(s):
 with (ROOT/'docs/TRANSLATIONAL_PROGRESS_LOG.md').open('a',encoding='utf-8') as f:f.write('\n'+s+'\n')
def write_json(obj,p,phase,desc):
 p=Path(p)
 if p.exists():raise FileExistsError(p)
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)),encoding='utf-8');register(p,phase,desc)
def write_csv(df,p,phase,desc):
 p=Path(p)
 if p.exists():raise FileExistsError(p)
 p.parent.mkdir(parents=True,exist_ok=True);df.to_csv(p,index=False);register(p,phase,desc)
def figure(fig,phase,name):
 d=ROOT/f'figures/translational/{phase}';d.mkdir(parents=True,exist_ok=True)
 stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');paths=[]
 for ext in ['png','pdf']:
  p=d/f'{name}_{stamp}.{ext}';fig.savefig(p,dpi=180,bbox_inches='tight');register(p,phase,name);paths.append(rel(p))
 return paths
def download(name,url,log):
 DATA.mkdir(parents=True,exist_ok=True);p=DATA/name;p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists():return p,{'url':url,'file':name,'status':'CACHED','bytes':p.stat().st_size,'sha256':sha(p)}
 try:
  with urllib.request.urlopen(url,timeout=120) as r:
   first=r.read(32768)
   if b'reCAPTCHA' in first or b'Checking your browser' in first:raise RuntimeError('CAPTCHA page; no bypass attempted')
   part=p.with_name(p.name+'.part')
   with part.open('wb') as o:
    o.write(first)
    while True:
     b=r.read(1<<20)
     if not b:break
     o.write(b)
  part.rename(p)
  rec={'url':url,'file':name,'status':'DOWNLOADED','bytes':p.stat().st_size,'sha256':sha(p)}
  log.info('download %s %s bytes',name,rec['bytes']);return p,rec
 except Exception as e:
  rec={'url':url,'file':name,'status':'BLOCKED','error':str(e),'manual_action':'Download this URL manually into F:/Thesis/data/tcga/'+name}
  log.warning('%s',rec);return None,rec
def bh(p):
 p=np.asarray(p,float);q=np.full(p.shape,np.nan);ok=np.isfinite(p);v=p[ok];n=len(v)
 if n:
  order=np.argsort(v);adj=np.minimum.accumulate((v[order]*n/np.arange(1,n+1))[::-1])[::-1];r=np.empty(n);r[order]=np.minimum(adj,1);q[ok]=r
 return q
