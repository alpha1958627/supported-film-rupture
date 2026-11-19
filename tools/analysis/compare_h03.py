"""Threshold crossings, registered TCW+mode and old saved bounds; default hmin=.3."""
from pathlib import Path
import argparse,json,hashlib,subprocess,base64,sys,os,csv,io
import numpy as np

parser=argparse.ArgumentParser();parser.add_argument('--workspace',required=True);parser.add_argument('--threshold',type=float,choices=(.3,.4),default=.3);args=parser.parse_args()
THRESHOLD=args.threshold;assert 0<THRESHOLD<.5
TAG='h%02d'%round(THRESHOLD*10);COUNT_KEY='n_reached_%02d'%round(THRESHOLD*10)
ROOT=Path(args.workspace).resolve();REPO=Path(__file__).resolve().parents[2];OUT=REPO/'dashboard'
snapshot=json.loads((ROOT/'snapshot.json').read_text())
items=[v for group in snapshot['groups'] for v in group.get('items',[]) if v['route']['kind']=='production']
assert len(items)==25
remote=r'''
import json,hashlib,io
from pathlib import Path
import numpy as np
items=ITEMS
threshold=THRESHOLD
be=[];plot_points={}
for v in items:
 p=Path(v['remote']);steps=v['accepted']
 with (p/'run/history.csv').open('rb') as f:raw=b''.join(f.readline() for _ in range(steps+1))
 assert hashlib.sha256(raw).hexdigest()==v['history_prefix_sha256']
 data=np.loadtxt(io.BytesIO(raw),delimiter=',',skiprows=1,usecols=(0,2),ndmin=2)
 assert len(data)==steps and data[-1,0]==v['time']
 hit=np.flatnonzero(data[:,1]<=threshold)
 record=dict(name=v['name'],height=v['height'],seed=v['seed'],reached=bool(len(hit)),observed_tau=v['time'],history_prefix_sha256=v['history_prefix_sha256'])
 if len(hit):
  i=int(hit[0]);before=data[i-1] if i else [0,v['height']]
  record.update(tau_lower=float(before[0]),tau_upper=float(data[i,0]),h_before=float(before[1]),h_after=float(data[i,1]),accepted_step=i+1)
  visible=data[:i+1]
 else:visible=data
 logt=np.log(visible[:,0]);span=max(float(logt[-1]-logt[0]),1e-30)
 bins=np.minimum(((logt-logt[0])/span*1200).astype(int),1199);keep=set([0,len(visible)-1])
 starts=np.r_[0,np.flatnonzero(np.diff(bins))+1];ends=np.r_[starts[1:],len(visible)]
 for lo,hi in zip(starts,ends):keep.update([int(lo),int(hi-1),int(lo+np.argmin(visible[lo:hi,1])),int(lo+np.argmax(visible[lo:hi,1]))])
 plot_points[v['name']]=visible[sorted(keep)].tolist()
 be.append(record)
old=[]
for task in range(20):
 p=Path('/public/home/ljh/basilisk_midpoint_formal_20261025/formal/task_%02d'%task)
 protocol=json.loads((p/'protocol.json').read_text());height=protocol['case']['height'];seed=protocol['case']['seed']
 assert protocol['cap']==.3 and protocol['phi_ref']==.01 and protocol['ls_ref']==.2
 raw=(p/'observations.jsonl').read_bytes();observations=[json.loads(line) for line in raw.splitlines()]
 hit=next(v for v in observations if height*v['hmin']<=threshold)
 old.append(dict(task=task,height=height,seed=seed,tau_first_saved_below=height*hit['time']/1e8,h_first_saved_below=height*hit['hmin'],observation_interval_native=300,observations_sha256=hashlib.sha256(raw).hexdigest(),interpretation='Upper bound on true first crossing; earlier unsaved downcrossings may exist. No interpolation or lower bracket is claimed.'))
print(json.dumps(dict(BE=be,Basilisk=old,plot_points=plot_points)))
'''
records=[]
for v in items:
 route=v['route'];d=v['latest'];records.append(dict(name=route['name'],remote=route['remote'],height=d['config']['initial'],seed=d['config']['seed'],accepted=d['accepted'],time=d['time'],history_prefix_sha256=v['history_prefix_sha256']))
assert len({v['route']['host'] for v in items})==1
code=remote.replace('ITEMS',repr(records)).replace('THRESHOLD',repr(THRESHOLD));command='echo '+base64.b64encode(code.encode()).decode()+' | base64 -d | '+items[0]['route']['python']+' -'
cachepath=ROOT/('crossings_'+TAG+'_raw.json');cached=json.loads(cachepath.read_text()) if cachepath.exists() else {}
oldrows={v['name']:v for v in cached.get('BE',[])}
same_snapshot='plot_points' in cached and len(oldrows)==25 and all(oldrows.get(v['name'],{}).get('history_prefix_sha256')==v['history_prefix_sha256'] and oldrows[v['name']]['observed_tau']==v['time'] for v in records)
if same_snapshot:data=cached
else:
 result=subprocess.run(['ssh','-o','BatchMode=yes',items[0]['route']['host'],command],capture_output=True,text=True,timeout=240)
 if result.returncode:raise RuntimeError(result.stderr)
 data=json.loads(result.stdout);cachepath.write_text(json.dumps(data,indent=2),encoding='utf-8')
sys.path.insert(0,str(REPO/'be'));from tcw_mode import TCWMode
models={h:TCWMode(h) for h in [.5,.6,.7,.8,.9]};summary=[]
for height,m in models.items():
 group=[v for v in data['BE'] if v['height']==height];be=[v['tau_upper'] for v in group if v['reached']];old=[v['tau_first_saved_below'] for v in data['Basilisk'] if v['height']==height]
 prediction=m.crossing(THRESHOLD);assert abs(m.predict(prediction)-THRESHOLD)<1e-8
 row=dict(height=height,n_observed=5,BE_mean_tau=float(np.mean(be)) if len(be)==5 else None,BE_SD=float(np.std(be,ddof=1)) if len(be)==5 else None,BE_event_tau=be,
  TCW_mode_tau=prediction,TCW_branch='TCW' if prediction<=m.switch else 'mode',Basilisk_n=len(old),Basilisk_mean_saved_upper_bound=float(np.mean(old)) if old else None,Basilisk_saved_upper_bounds=old)
 row[COUNT_KEY]=len(be);summary.append(row)
report=dict(snapshot_utc=snapshot['updated_utc'],threshold=THRESHOLD,phi=.01,C_rms=3.5,summary=summary,events={k:data[k] for k in ['BE','Basilisk']},
 note='Same previously verified 25-path snapshot, reclassified at first hmin/h0<=%g. BE uses all accepted-step records. Old Basilisk observation spacing is 300 native-time units, so blue markers are upper bounds, not exact first-passage times. Simulation jobs retain their original .005 stop setting; this is an analysis criterion.'%THRESHOLD,
 theory='One fixed no-slip white6 continuum TCW+mode construction, fixed 256² sum. No fitted C or seed-dependent theory. This is an extremal-envelope prediction, not a derived stochastic MFPT.',
 normalization='New BE directly outputs tau; old tau=height*t_native/1e8, h/h0=height*H_native. Old slip=.2, new no slip. Exact old batch has no height=.5.')
(OUT/(TAG+'_comparison.json')).write_text(json.dumps(report,indent=2),encoding='utf-8')
for name,rows in [('BE_'+TAG+'_events',data['BE']),('Basilisk_'+TAG+'_saved_bounds',data['Basilisk'])]:
 stream=io.StringIO();writer=csv.DictWriter(stream,fieldnames=list(dict.fromkeys(k for row in rows for k in row)),lineterminator='\n');writer.writeheader();writer.writerows(rows);(OUT/(name+'.csv')).write_bytes(stream.getvalue().encode('utf-8-sig'))
os.environ['MPLCONFIGDIR']=str(ROOT/'mpl_cache');sys.path.insert(0,str(REPO/'postprocess'))
from academic_style import apply_academic_style
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
apply_academic_style();fig=plt.figure(figsize=(8,6));ax=fig.add_axes([.16,.17,.81,.55])
height_grid=np.linspace(.5,.9,17);theory=[TCWMode(h).crossing(THRESHOLD) for h in height_grid]
ax.plot(height_grid,theory,color='k',lw=2)
for height in models:
 group=[v for v in data['BE'] if v['height']==height];hits=[v for v in group if v['reached']]
 for i,v in enumerate(hits):ax.plot(height-.012+(i-2)*.002,v['tau_upper'],'o',mfc='white',mec='#D62728',ms=5)
 if len(hits)==5:
  values=np.array([v['tau_upper'] for v in hits]);ax.errorbar(height-.012,values.mean(),yerr=values.std(ddof=1),color='#D62728',fmt='s',ms=5,capsize=4,lw=1.3)
 for i,v in enumerate([v for v in group if not v['reached']]):
  x=height-.012+(i-2)*.002;y=v['observed_tau'];ax.plot(x,y,'^',mfc='white',mec='#D62728',ms=5);ax.annotate('',xy=(x,y*1.5),xytext=(x,y),arrowprops=dict(arrowstyle='->',color='#D62728',lw=.7))
 old=[v for v in data['Basilisk'] if v['height']==height]
 for i,v in enumerate(old):
  x=height+.012+(i-2)*.002;y=v['tau_first_saved_below'];ax.plot(x,y,'v',mfc='white',mec='#1F77B4',ms=5);ax.annotate('',xy=(x,y/1.4),xytext=(x,y),arrowprops=dict(arrowstyle='->',color='#1F77B4',lw=.7))
ax.set(xlabel=r'$h_i/h_0$',ylabel=r'$\tau_{%g}$'%THRESHOLD,yscale='log',xlim=(.46,.94));ax.set_xticks([.5,.6,.7,.8,.9]);ax.grid(False,which='both')
handles=[Line2D([],[],color='k',lw=2,label='TCW + mode'),Line2D([],[],color='#D62728',marker='o',mfc='white',ls='none',label='BE: first crossing'),Line2D([],[],color='#1F77B4',marker='v',mfc='white',ls='none',label='Basilisk .3: upper bounds'),Line2D([],[],color='#D62728',marker='^',mfc='white',ls='none',label='BE: incomplete')]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.995),ncol=1,labelspacing=.15)
stem=OUT/('TCW_mode_vs_'+TAG);fig.savefig(str(stem)+'.pdf');fig.savefig(str(stem)+'.png',dpi=600);fig.savefig(str(stem)+'_preview.png',dpi=120);plt.close(fig)
fig=plt.figure(figsize=(8,6));ax=fig.add_axes([.16,.17,.81,.55]);lookup={v['name']:v for v in data['BE']}
colors=['#000000','#D62728','#1F77B4','#9467BD','#2CA02C'];markers=['o','s','^','D','v']
for (height,m),color in zip(models.items(),colors):
 for item in [v for v in items if v['latest']['config']['initial']==height]:
  event=lookup[item['route']['name']];ar=np.array(data['plot_points'][item['route']['name']])
  if event['reached']:
   ar=ar[ar[:,0]<event['tau_lower']]
   ar=np.vstack((ar,[event['tau_lower'],event['h_before']],[event['tau_upper'],event['h_after']]))
  ar=ar[ar[:,0]>0];ax.plot(ar[:,0],ar[:,1],color=color,lw=.6,alpha=.8)
  ax.plot(ar[-1,0],ar[-1,1],'o' if event['reached'] else 's',mfc='white',mec=color,ms=4)
 end=m.crossing(THRESHOLD);times=np.geomspace(1e-12,end,500);prediction=[m.predict(t) for t in times]
 ax.plot(times,prediction,color=color,lw=1.8,ls='--')
ax.axhline(THRESHOLD,color='.5',ls=':',lw=.8);ax.set(xlabel=r'$\tau$',ylabel=r'$h_{\min}/h_0$',xscale='log',xlim=(1e-10,max(m.crossing(THRESHOLD) for m in models.values())*1.1),ylim=(THRESHOLD-.03,.94));ax.grid(False,which='both')
handles=[Line2D([],[],color=c,marker=s,mfc='white',lw=1,label=r'$h_i/h_0=%.1f$'%h) for h,c,s in zip(models,colors,markers)]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.995),ncol=3,columnspacing=.6,handlelength=1.2)
fig.text(.16,.77,r'Solid: BE; dashed: TCW + mode',va='bottom')
stem=OUT/('hmin_first_'+TAG);fig.savefig(str(stem)+'.pdf');fig.savefig(str(stem)+'.png',dpi=600);fig.savefig(str(stem)+'_preview.png',dpi=120);plt.close(fig)
lines=['# 按首次到hmin/h0=%g定义事件：TCW+mode vs BE / Basilisk .3'%THRESHOLD,'', '同一核验快照UTC：'+snapshot['updated_utc']+'。本页重分类统计阈值，原运行端仍按.005保存后续轨迹。','',
 '|初厚|BE到%g数量|BE平均tau（5条全到才给均值）|TCW+mode到%g|Basilisk首次保存低于%g的时间上界均值|'%(THRESHOLD,THRESHOLD,THRESHOLD),'|---|---:|---:|---:|---:|']
f=lambda x:'—' if x is None else '%.8g'%x
for row in summary:lines.append('|%.1f|%d/5|%s|%s|%s|'%(row['height'],row[COUNT_KEY],f(row['BE_mean_tau']),f(row['TCW_mode_tau']),f(row['Basilisk_mean_saved_upper_bound'])))
early='、'.join('%.1f'%row['height'] for row in summary if row['TCW_branch']=='TCW') or '无'
lines+=['', '黑线：固定C_rms=3.5、phi=.01、无滑移white6 TCW+mode；不是对随机首达平均的推导。TCW分支先到阈值的厚度：'+early+'；这些组直接求TCW分支交点，其他组用mode分支。', '',
 '红色圆点为每条BE首次接受步下穿，红方块/误差棒为5条均值±样本SD；红上箭头为尚未到阈值的观察下界。蓝色倒三角/下箭头为旧Basilisk首次保存低于阈值的时间上界：旧观测每300原生时间输出，可能遗漏更早的下穿与回升，不把相邻观测当严格首达括区。', '',
 '旧tau=height*t_native/1e8；新tau直接输出。旧有滑移.2、新无滑移，不能当同条件数值误差验证。旧精确批次没有.5组。', '',
 '[BE逐条事件](BE_'+TAG+'_events.csv) · [旧保存上界](Basilisk_'+TAG+'_saved_bounds.csv) · [完整参数与来源SHA]('+TAG+'_comparison.json)']
(OUT/(TAG+'_comparison.md')).write_text('\n'.join(lines),encoding='utf-8')
page=OUT/'index.html';html=page.read_text(encoding='utf-8')
total=sum(v[COUNT_KEY] for v in summary)
section='<section id="%s"><h2>首次到hmin/h0=%g：TCW+mode vs 数值轨迹</h2><p class="note">按此判据，%d/25条已到达，%d条尚未到达。旧Basilisk仅能由稀疏保存数据得到首达上界，原0.005分析保留。</p>'%(TAG,THRESHOLD,total,25-total)
section+='<a href="TCW_mode_vs_'+TAG+'.png"><img style="max-width:960px" src="TCW_mode_vs_'+TAG+'_preview.png"></a><p><a href="'+TAG+'_comparison.md">统计与定义</a> · <a href="BE_'+TAG+'_events.csv">逐条BE事件</a></p><a href="hmin_first_'+TAG+'.png"><img style="max-width:960px" src="hmin_first_'+TAG+'_preview.png"></a></section>'
if 'id="'+TAG+'"' in html:
 start=html.index('<section id="'+TAG+'">');end=html.index('</section>',start)+len('</section>');html=html[:start]+html[end:]
page.write_text(html.replace('</main>',section+'</main>'),encoding='utf-8')
print(json.dumps(summary,indent=2))
