"""Read selected no-slip MPI jobs and render the source-backed hmin board."""
from pathlib import Path
import argparse,json,subprocess,shlex,hashlib,sys,os,csv,io
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--workspace',required=True);args=p.parse_args()
ROOT=Path(args.workspace);REPO=Path(__file__).resolve().parent;OUT=REPO/'dashboard/basilisk_noslip';OUT.mkdir(exist_ok=True)
campaign=json.loads((ROOT/'campaign.json').read_text());selected=json.loads((ROOT/'selected.json').read_text());records=[]
code=r'''
from pathlib import Path
import json,hashlib,io
import numpy as np
p=Path.cwd();row={}
for name in ('validation.json','validation_failed.json','structural_pass.json'):
 if (p/name).exists():row[name]=json.loads((p/name).read_text())
if (p/'run/latest.json').exists():
 d=json.loads((p/'run/latest.json').read_text());row['latest']=d
 checkpoint=p/'run/checkpoint.npz';assert hashlib.sha256(checkpoint.read_bytes()).hexdigest()==d['checkpoint_sha256']
 with (p/'run/history.csv').open('rb') as f:raw=b''.join(f.readline() for _ in range(d['accepted']+1))
 assert hashlib.sha256(raw).hexdigest()==d['history_sha256']
 if d['accepted']:
  a=np.loadtxt(io.BytesIO(raw),delimiter=',',skiprows=1,usecols=(0,2),ndmin=2);assert len(a)==d['accepted'] and a[-1,0]==d['time']
  # Preserve extrema on a log-time sampling mesh, including earliest transients.
  keep={0,len(a)-1};bins=np.searchsorted(np.geomspace(max(a[0,0],1e-20),max(a[-1,0],a[0,0])*1.00001,401),a[:,0])
  for b in np.unique(bins):
   ix=np.flatnonzero(bins==b);keep.update((int(ix[0]),int(ix[-1]),int(ix[np.argmin(a[ix,1])]),int(ix[np.argmax(a[ix,1])])))
  row['curve']=a[sorted(keep)].tolist()
 else:row['curve']=[]
print(json.dumps(row))
'''
routes=json.loads((ROOT/'candidates.json').read_text())
pending={r['name']:dict(winner=r) for r in reversed(routes) if r['name'] not in selected}
for name,choice in dict(selected,**pending).items():
 route=choice['winner'];cmd='cd '+shlex.quote(route['remote'])+' && PYTHONPATH='+shlex.quote(route.get('pythonpath',''))+' '+route['python']+' -c '+shlex.quote(code)
 r=subprocess.run(['ssh',route['host'],'bash -lc '+shlex.quote(cmd)],capture_output=True,text=True,timeout=90)
 item=dict(name=name,server=route['server'],job=route['job'],kind=route['kind'],selected=name in selected)
 if r.returncode:item['read_error']=r.stderr[-500:]
 else:item.update(json.loads(r.stdout))
 # Independent scheduler state: a threshold event and Slurm completion differ.
 cmd='squeue -h -j '+route['job']+' -o "%T|%S|%N"; sacct -nP -j '+route['job']+' -X -o State,Start,NodeList'
 r=subprocess.run(['ssh',route['host'],'bash -lc '+shlex.quote(cmd)],capture_output=True,text=True,timeout=30)
 item['scheduler']=r.stdout.strip().splitlines()[0] if r.returncode==0 and r.stdout.strip() else 'QUERY_ERROR'
 records.append(item)
snapshot=dict(utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),records=records)
(ROOT/'bm_snapshot.json').write_text(json.dumps(snapshot,indent=2))
# Publish metadata and sampled hmin only; no raw fields, checkpoint or full history.
public=[]
for r in records:
 row=dict(r)
 if 'latest' in row:
  row['latest']={k:v for k,v in row['latest'].items() if k not in ('config','ledger')}
 public.append(row)
(OUT/'data.json').write_text(json.dumps(dict(utc=snapshot['utc'],records=public),indent=2))
sys.path.insert(0,str(REPO/'basilisk_midpoint'));from tcw_mode import TCWMode
os.environ['MPLCONFIGDIR']=str(ROOT/'mpl_cache');sys.path.insert(0,str(REPO/'postprocess'))
from academic_style import apply_academic_style
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
apply_academic_style();heights=[.5,.6,.7,.8,.9];colors=['k','#D62728','#1F77B4','#9467bd','#2ca02c'];markers=['o','s','^','D','v'];models={h:TCWMode(h) for h in heights}
def save(fig,name):
 fig.savefig(OUT/(name+'.pdf'));fig.savefig(OUT/(name+'.png'),dpi=600);fig.savefig(OUT/(name+'_preview.png'),dpi=120);plt.close(fig)
def figure():
 fig=plt.figure(figsize=(8,6));ax=fig.add_axes([.16,.21,.81,.49]);ax.set(xscale='log',xlabel=r'$\tau$',ylabel=r'$h_{\min}/h_0$',ylim=(.37,.96));ax.grid(False,which='both');ax.axhline(.4,color='0.5',ls=':',lw=1)
 fig.legend(handles=[Line2D([],[],color=c,marker=m,mfc='white',label=r'$h_i/h_0=%.1f$'%h) for h,c,m in zip(heights,colors,markers)],loc='upper center',bbox_to_anchor=(.55,1.01),ncol=3,columnspacing=.7,handlelength=1.3,labelspacing=.1)
 return fig,ax
fig,ax=figure();end=max(m.crossing(.4) for m in models.values());ax.set_xlim(1e-10,end*1.1)
production=[r for r in records if r['kind']=='production'];summary=[]
for h,c in zip(heights,colors):
 m=models[h];t=np.geomspace(1e-10,m.crossing(.4),450);ax.plot(t,[m.predict(v) for v in t],ls='--',color=c,lw=1.7)
 group=[r for r in production if r.get('latest',{}).get('hmin') is not None and next(v for v in campaign['cases'] if v['name']==r['name'])['initial']==h]
 for r in group:
  a=np.array(r.get('curve',[]))
  if len(a):ax.plot(a[:,0],a[:,1],color=c,lw=.7,alpha=.65)
  d=r['latest'];ax.plot(max(d['time'],1e-10),d['hmin'],marker='o' if d['status']=='THRESHOLD' else 's',mfc='white',mec=c,ms=4)
 summary.append(dict(initial=h,production_submitted=sum(next(v for v in campaign['cases'] if v['name']==r['name'])['initial']==h for r in production),threshold=sum(r['latest']['status']=='THRESHOLD' for r in group),theory_tau=m.crossing(.4),minimum=min((r['latest']['hmin'] for r in group),default=None)))
fig.text(.55,.765,r'Basilisk: solid; TCW + mode: dashed',ha='center');fig.text(.55,.025,r'No slip, $\varphi=0.01$, $C_{\rm rms}=5$',ha='center');save(fig,'hmin_overview')
lines=['# Basilisk 无滑移 hmin 看板','', '更新时间 UTC '+snapshot['utc'],'','64²、每条4 MPI核；φ=.01、m=u³、A=13.159、C_rms=5。实线仅为已保存的新模拟；虚线为64²连续Fourier谱TCW+mode，无白噪声插值滤波。固定绝对阈值u=h/h0=.4。尚未计算的区间仅有模型虚线。','', '|工况|服务器|作业|类型|调度状态|u_min|tau|','|---|---|---|---|---|---:|---:|']
for r in records:
 d=r.get('latest',{});lines.append('|%s|%s|%s|%s|%s|%s|%s|'%(r['name'],r['server'],r['job'],r['kind'],r['scheduler'],d.get('hmin','—'),d.get('time','—')))
lines+=['','![hmin](hmin_overview.png)','','短测逐厚度通过后启动其5条路径；未通过不启动。6小时资源窗保存检查点续排，不代表物理完成。原BE任务不变。源码见../../basilisk_midpoint。']
(OUT/'README.md').write_text('\n'.join(lines),encoding='utf-8');(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
(OUT/'index.html').write_text('<meta charset="utf-8"><title>Basilisk 无滑移看板</title><h1>Basilisk 无滑移：25条轨迹</h1><p>UTC '+snapshot['utc']+'</p><img style="max-width:100%" src="hmin_overview_preview.png"><p><a href="README.md">任务表、判据与状态</a> · <a href="data.json">可核验数据</a></p>',encoding='utf-8')
print(json.dumps(summary,indent=2))
