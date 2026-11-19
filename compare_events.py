"""Compare saved first passages in common units and at the same height threshold."""
from pathlib import Path
import argparse,json,csv,hashlib,subprocess,base64,io,sys,os,datetime
import numpy as np

parser=argparse.ArgumentParser();parser.add_argument('--workspace',required=True);args=parser.parse_args()
ROOT=Path(args.workspace).resolve();REPO=Path(__file__).resolve().parent;OUT=REPO/'dashboard'
snapshot=json.loads((ROOT/'snapshot.json').read_text())
items=[v for g in snapshot['groups'] for v in g.get('items',[]) if v['route']['kind']=='production']
assert len(items)==25 and not any(v.get('numerical_gate_violation') for v in items)
completed=[v for v in items if v['latest']['status']=='THRESHOLD']
for v in completed:
 d=v['latest'];e=d['event'];assert e['hmin_before']>.005 and 0<e['hmin_after']<=.005 and e['upper']==d['time'] and e['upper']-e['lower']<=1e-6/805306368
cachepath=ROOT/'matched_events.json';cache=json.loads(cachepath.read_text()) if cachepath.exists() else {}
needed=[v for v in completed if cache.get(v['route']['name'],{}).get('checkpoint_sha256')!=v['latest']['checkpoint_sha256']]
remote=r'''
from pathlib import Path
import json,hashlib,io
import numpy as np
items=ITEMS
result={}
for item in items:
 p=Path(item['route']['remote']);d=item['latest'];threshold=.01*d['config']['initial']
 assert hashlib.sha256((p/'run/checkpoint.npz').read_bytes()).hexdigest()==d['checkpoint_sha256']
 with (p/'run/history.csv').open('rb') as f:raw=b''.join(f.readline() for _ in range(d['accepted']+1))
 assert hashlib.sha256(raw).hexdigest()==item['history_prefix_sha256']
 data=np.loadtxt(io.BytesIO(raw),delimiter=',',skiprows=1,ndmin=2)
 assert len(data)==d['accepted'] and data[-1,0]==d['time']
 hits=np.flatnonzero(data[:,2]<=threshold);assert len(hits)
 j=int(hits[0]);before=data[j-1] if j else [0,0,d['config']['initial']]
 result[item['route']['name']]=dict(checkpoint_sha256=d['checkpoint_sha256'],history_prefix_sha256=item['history_prefix_sha256'],threshold=threshold,tau_lower=float(before[0]),tau_upper=float(data[j,0]),h_before=float(before[2]),h_after=float(data[j,2]),accepted_step=j+1)
print(json.dumps(result))
'''
for host in sorted({v['route']['host'] for v in needed}):
 own=[v for v in needed if v['route']['host']==host]
 compact=[dict(route={k:v['route'][k] for k in ['remote','name']},latest={k:v['latest'][k] for k in ['accepted','time','checkpoint_sha256','config']},history_prefix_sha256=v['history_prefix_sha256']) for v in own]
 code=remote.replace('ITEMS',repr(compact));command='echo '+base64.b64encode(code.encode()).decode()+' | base64 -d | '+own[0]['route']['python']+' -'
 r=subprocess.run(['ssh','-o','BatchMode=yes',host,command],capture_output=True,text=True,timeout=180)
 if r.returncode:raise RuntimeError(r.stderr)
 cache.update(json.loads(r.stdout))
cachepath.write_text(json.dumps(cache,indent=2),encoding='utf-8')

source=REPO/'reference/basilisk_cap03.csv';provenance=json.loads((REPO/'reference/provenance.json').read_text())
assert hashlib.sha256(source.read_bytes()).hexdigest()==provenance['derived_reference_sha256']
old=list(csv.DictReader(source.open(encoding='utf-8-sig')))
for v in old:
 assert v['outcome']=='RUPTURE_BRACKET' and v['raw_verified']=='True'
 assert np.isclose(float(v['tau']),float(v['height'])*float(v['native_time'])/1e8,rtol=1e-14)
accounting={}
for group in snapshot['groups']:
 for line in group['scheduler'].splitlines():
  a=line.split('|')
  if len(a)>=7:accounting[a[0]]=a
def seconds(value):
 days,clock=value.split('-',1) if '-' in value else ('0',value)
 h,m,s=clock.split(':');return int(days)*86400+int(h)*3600+int(m)*60+float(s)

events=[]
for v in completed:
 d=v['latest'];route=v['route'];match=cache[route['name']];ac=accounting[route['job']]
 assert ac[1]=='COMPLETED'
 events.append(dict(name=route['name'],height=d['config']['initial'],replicate=d['config']['replicate'],seed=d['config']['seed'],job=route['job'],server=route['server'],physical_cores=4,
  tau_at_005=d['event']['upper'],stop_lower=d['event']['lower'],stop_upper=d['event']['upper'],stop_h_before=d['event']['hmin_before'],stop_h_after=d['event']['hmin_after'],
  matched_threshold=match['threshold'],matched_tau_lower=match['tau_lower'],matched_tau_upper=match['tau_upper'],matched_h_before=match['h_before'],matched_h_after=match['h_after'],
  slurm_wall_seconds=seconds(ac[3]),main_wall_seconds=d['actual_wall_seconds'],actual_start_beijing=ac[5],actual_end_beijing=ac[6],checkpoint_sha256=d['checkpoint_sha256'],history_prefix_sha256=v['history_prefix_sha256']))
events.sort(key=lambda v:(v['height'],v['replicate']))
def save_csv(path,rows):
 stream=io.StringIO();w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows);path.write_bytes(stream.getvalue().encode('utf-8-sig'))
if events:save_csv(OUT/'completed_events.csv',events)
summary=[]
for height in (.5,.6,.7,.8,.9):
 group=[v for v in items if v['latest']['config']['initial']==height];be=[v for v in events if v['height']==height];ref=[v for v in old if float(v['height'])==height]
 row=dict(height=height,n_completed=len(be),n_running=sum(v['latest']['status']=='RUNNING' for v in group),n_numerical_failure=sum(v['latest']['status']=='NUMERICAL_FAILURE' for v in group),n_basilisk=len(ref))
 for field,key in [('tau_at_005','BE_tau_005'),('matched_tau_upper','BE_tau_matched'),('slurm_wall_seconds','BE_wall_seconds')]:
  values=[v[field] for v in be]
  row[key+'_mean']=float(np.mean(values)) if len(be)==5 else None
  row[key+'_SD']=float(np.std(values,ddof=1)) if len(be)==5 else None
  row[key+'_range']=[min(values),max(values)] if values else None
 row['Basilisk_tau_mean']=float(np.mean([float(v['tau']) for v in ref])) if len(ref)==5 else None
 row['Basilisk_tau_SD']=float(np.std([float(v['tau']) for v in ref],ddof=1)) if len(ref)==5 else None
 row['Basilisk_wall_seconds_mean']=float(np.mean([float(v['slurm_wall_seconds']) for v in ref])) if len(ref)==5 else None
 row['BE_over_Basilisk_tau_ratio']=row['BE_tau_matched_mean']/row['Basilisk_tau_mean'] if len(be)==len(ref)==5 else None
 row['BE_over_Basilisk_wall_ratio']=row['BE_wall_seconds_mean']/row['Basilisk_wall_seconds_mean'] if len(be)==len(ref)==5 else None
 alive=[v['latest'] for v in group if v['latest']['status']=='RUNNING']
 row['running_tau_range']=[min(v['time'] for v in alive),max(v['time'] for v in alive)] if alive else None
 row['running_hmin_range']=[min(v['hmin'] for v in alive),max(v['hmin'] for v in alive)] if alive else None
 summary.append(row)
report=dict(snapshot_utc=snapshot['updated_utc'],groups=summary,completed_events=events,reference=provenance,
 normalization=dict(BE='native tau output; h_i/h_0=.5-.9',Basilisk='tau=(h_i/h_0)*t_native/1e8; h/h_0=(h_i/h_0)*H_native'),
 comparison_threshold='Matched at hmin=.01*h_initial. BE production stop remains .005.',
 cautions=['No h=.5 data exists in the exact Basilisk .3 reference batch.', 'BE no-slip versus Basilisk slip=.2; BE versus midpoint, 4 physical OpenMP cores versus 16 MPI ranks; not a matched-method benchmark.', 'Error bars are sample SD, n=5. Running paths are right-censored, not rupture events or event-subset means.'])
(OUT/'event_comparison.json').write_text(json.dumps(report,indent=2),encoding='utf-8')

os.environ['MPLCONFIGDIR']=str(ROOT/'mpl_cache');sys.path.insert(0,str(REPO/'postprocess'))
from academic_style import apply_academic_style
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
apply_academic_style()
fig=plt.figure(figsize=(8,6));ax=fig.add_axes([.16,.17,.81,.57])
for height in (.5,.6,.7,.8,.9):
 be=[v for v in events if v['height']==height];ref=[v for v in old if float(v['height'])==height]
 for rows,offset,color,key in [(be,-.012,'k','matched_tau_upper'),(ref,.012,'#D62728','tau')]:
  values=np.array([float(v[key]) for v in rows])
  if len(values):ax.scatter(height+offset+np.linspace(-.005,.005,len(values)),values,s=27,facecolors='none',edgecolors=color,zorder=3)
  if len(values)==5:ax.errorbar(height+offset,values.mean(),yerr=values.std(ddof=1),color=color,fmt='s',ms=5,capsize=4,lw=1.6,zorder=4)
 for i,v in enumerate([v for v in items if v['latest']['config']['initial']==height and v['latest']['status']=='RUNNING']):
  d=v['latest'];assert v['min_hmin']>.01*height
  x=height-.012+(i-2)*.003;t=d['time'];ax.plot(x,t,'^',mfc='white',mec='k',ms=5)
  ax.annotate('',xy=(x,t*1.35),xytext=(x,t),arrowprops=dict(arrowstyle='->',color='k',lw=.8))
ax.set(xlabel=r'$h_i/h_0$',ylabel=r'First-passage time $\tau$',yscale='log',xlim=(.46,.94));ax.set_xticks([.5,.6,.7,.8,.9]);ax.grid(False,which='both')
handles=[Line2D([],[],color='k',marker='o',mfc='white',ls='none',label='BE, no slip'),Line2D([],[],color='#D62728',marker='o',mfc='white',ls='none',label=r'Basilisk .3, $\ell=0.2$'),Line2D([],[],color='k',marker='^',mfc='white',ls='none',label='BE: incomplete')]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.995),ncol=1,labelspacing=.25)
stem=OUT/'rupture_comparison'
fig.savefig(str(stem)+'.pdf');fig.savefig(str(stem)+'.png',dpi=600);fig.savefig(str(stem)+'_preview.png',dpi=120);plt.close(fig)
lines=['# BE与Basilisk .3：最新破碎时间比较','', '采集快照UTC：'+snapshot['updated_utc']+'。下表的完成指BE实测首次跨hmin=.005，未到阈值不计完成。','',
 '|初厚|BE完成/5|运行|BE到.005平均tau|BE到共同阈值平均tau|Basilisk .3平均tau|共同阈值时间比BE/旧|','|---|---:|---:|---:|---:|---:|---:|']
def f(v):return '—' if v is None else '%.8g'%v
for v in summary:lines.append('|%.1f|%d|%d|%s|%s|%s|%s|'%(v['height'],v['n_completed'],v['n_running'],f(v['BE_tau_005_mean']),f(v['BE_tau_matched_mean']),f(v['Basilisk_tau_mean']),f(v['BE_over_Basilisk_tau_ratio'])))
lines+=['', '共同阈值为初始厚度的1%：例如.6组为hmin=.006。BE原停止阈值仍.005；共同阈值从完整、检查点核验的历史中提取首次下穿及前后接受步括区。旧tau=(h_i/h_0)*t_native/1e8；新tau直接输出。原Basilisk .3数据没有.5组。','',
 '图中细空心圆是每条事件，方块及误差棒为5条均值±样本SD；空心三角和向上箭头是未完成BE轨迹的观察下界。未完成组不计算事件均值。', '',
 '新BE无滑移、4物理OpenMP核；旧Basilisk滑移.2、随机中点、16MPI，空间吸引项离散也不同。时间比不是同工况算法加速比。', '',
 '|BE轨迹|到.005的tau|到共同阈值的tau|实际墙钟分钟|','|---|---:|---:|---:|']
for v in events:lines.append('|%s|%.10g|%.10g|%.3f|'%(v['name'],v['tau_at_005'],v['matched_tau_upper'],v['slurm_wall_seconds']/60))
(OUT/'event_comparison.md').write_text('\n'.join(lines),encoding='utf-8')
page=OUT/'index.html';text=page.read_text(encoding='utf-8')
text=text.replace('<!-- event-comparison -->','')
if '<section id="event-comparison">' in text:text=text.split('<section id="event-comparison">')[0]+'</main></html>'
addition='<section id="event-comparison"><h2>已完成轨迹与Basilisk .3：共同无量纲时间、共同阈值</h2><p class="note">旧组有滑移.2，本组无滑移；误差棒为5条样本SD，箭头为未完成下界。0.5没有对应旧组。</p><a href="rupture_comparison.png"><img style="max-width:960px" src="rupture_comparison_preview.png"></a><p><a href="event_comparison.md">统计与说明</a> · <a href="completed_events.csv">逐条事件CSV</a> · <a href="event_comparison.json">完整统计和来源</a></p></section>'
page.write_text(text.replace('</main>',addition+'</main>'),encoding='utf-8')
print(json.dumps(summary,indent=2))
