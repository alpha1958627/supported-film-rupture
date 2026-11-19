"""All-seed hmin figures and a source-backed, portable HTML status board."""
from pathlib import Path
import json,sys,argparse,datetime,html,csv,io,os
import numpy as np

p=argparse.ArgumentParser();p.add_argument('--workspace',required=True);args=p.parse_args()
ROOT=Path(args.workspace).resolve();REPO=Path(__file__).resolve().parents[2];OUT=REPO/'dashboard';OUT.mkdir(exist_ok=True)
os.environ['MPLCONFIGDIR']=str(ROOT/'mpl_cache')
sys.path.insert(0,str(REPO/'be'));from tcw_mode import TCWMode
sys.path.insert(0,str(REPO/'postprocess'))
from academic_style import apply_academic_style
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
apply_academic_style()
snapshot=json.loads((ROOT/'snapshot.json').read_text());cases=snapshot['campaign']['cases']
items={v['route']['name']:v for group in snapshot['groups'] for v in group.get('items',[])}
divisors={item['validation.json']['initial']:item['validation.json']['selected_divisor'] for item in items.values() if item.get('validation.json',{}).get('passed')}
states={};starts={};ends={}
for group in snapshot['groups']:
 for line in group.get('scheduler','').splitlines():
  fields=line.split('|')
  if len(fields)>=4:states.setdefault(fields[0],fields)
  if len(fields)>=7:starts[fields[0]]=fields[5];ends[fields[0]]=fields[6]
colors=['#000000','#D62728','#1F77B4','#9467BD','#2CA02C'];markers=['o','s','^','D','v'];heights=[.5,.6,.7,.8,.9]
models={h:TCWMode(h) for h in heights};theory={};records=[]
pilot_tau=.0002370080211070061;pilot_hours=(2*3600+54*60+2)/3600
for height,model in models.items():
 end=model.crossing(.005);times=np.unique(np.r_[0,np.geomspace(1e-14,end,180),np.linspace(0,end,400),model.switch])
 prediction=np.array([model.predict(v) for v in times]);theory[height]=np.c_[times,prediction]
 np.savetxt(OUT/('tcw_h%02d.csv'%round(height*10)),theory[height],delimiter=',',header='tau,hmin',comments='')
for cfg in cases:
 item=items.get(cfg['name'],{});route=item.get('route',{});latest=item.get('latest',{});job=route.get('job','');state=states.get(job,['','NOT_SUBMITTED','',''])
 status='WAITING_VALIDATION' if not job else state[1]
 if latest.get('status')=='THRESHOLD':status='THRESHOLD'
 if latest.get('status')=='NUMERICAL_FAILURE' or item.get('numerical_gate_violation'):status='NUMERICAL_FAILURE'
 divisor=latest.get('config',{}).get('dt_divisor',divisors.get(cfg['initial']))
 model_hours=pilot_hours*models[cfg['initial']].crossing(.005)/models[.7].crossing(.005)*(divisor or 8)/8
 row=dict(name=cfg['name'],initial=cfg['initial'],replicate=cfg['replicate'],seed=cfg['seed'],server=route.get('server',''),job=job,node=state[2],status=status,
  tau=latest.get('time'),hmin=latest.get('hmin'),minimum_recorded=item.get('min_hmin'),accepted=latest.get('accepted'),rejected=latest.get('rejected'),
  actual_start_beijing=starts.get(job),actual_end_beijing=ends.get(job),physical_cores=4,source_sha256=latest.get('source_sha256'),
  mass=latest.get('relative_mass'),residual=item.get('max_full_residual'),checkpoint_utc=latest.get('updated_utc'),checkpoint_sha256=latest.get('checkpoint_sha256'),history_prefix_sha256=item.get('history_prefix_sha256'),
  dt_divisor=divisor,actual_wall_hours=latest.get('actual_wall_seconds',0)/3600 if latest else None,
  planning_hours=model_hours,planning_basis='TCW .005 time ratio scaled by independent h=.7 pilot 2h54m and selected dt divisor/8; not an FPT estimate or confidence interval')
 event=latest.get('event') or {};row.update(event_lower=event.get('lower'),event_upper=event.get('upper'),event_h_before=event.get('hmin_before'),event_h_after=event.get('hmin_after'))
 if starts.get(job) and starts[job] not in ('Unknown','None'):
  beginning=datetime.datetime.fromisoformat(starts[job]).replace(tzinfo=datetime.timezone(datetime.timedelta(hours=8)))
  row['planning_end_beijing']=(beginning+datetime.timedelta(hours=model_hours)).isoformat()
 if status=='THRESHOLD':row['TCW_work_remaining_hours']=0
 elif status!='NUMERICAL_FAILURE' and row['tau'] and row['actual_wall_hours']:
  row['measured_tau_per_actual_hour']=row['tau']/row['actual_wall_hours']
  row['TCW_work_remaining_hours']=max(0,models[cfg['initial']].crossing(.005)-row['tau'])/row['measured_tau_per_actual_hour']
 records.append(row)
csvfile=io.StringIO();writer=csv.DictWriter(csvfile,fieldnames=list(dict.fromkeys(k for row in records for k in row)),lineterminator='\n');writer.writeheader();writer.writerows(records);(OUT/'status.csv').write_bytes(csvfile.getvalue().encode('utf-8-sig'))

def axes():
 fig=plt.figure(figsize=(8,6));ax=fig.add_axes([.16,.17,.81,.56]);ax.set(xlabel=r'$\tau$',ylabel=r'$h_{\min}/h_0$',ylim=(0,.94));ax.grid(False,which='both');ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0));return fig,ax
def save(fig,name):
 fig.savefig(OUT/(name+'.pdf'));fig.savefig(OUT/(name+'.png'),dpi=600);fig.savefig(OUT/(name+'_preview.png'),dpi=120);plt.close(fig)
def curve(ax,cfg,color,style='-',lw=.6):
 points=items.get(cfg['name'],{}).get('points',[])
 if points:
  ar=np.vstack(([0,cfg['initial']],points));ax.plot(ar[:,0],ar[:,1],color=color,ls=style,lw=lw,alpha=.8)
  last=items[cfg['name']]['latest'];symbol='o' if last['status']=='THRESHOLD' else 'x' if last['status']=='NUMERICAL_FAILURE' else 's'
  ax.plot(ar[-1,0],ar[-1,1],symbol,color=color,mfc='white',ms=4)
  return ar[-1,0]
 return 0.
fig,ax=axes();maxseen=0.
for height,color,marker in zip(heights,colors,markers):
 for cfg in cases:
  if cfg['initial']==height:maxseen=max(maxseen,curve(ax,cfg,color))
 # Full model extent is available in a separate forecast figure; overview uses observed window.
 plotend=maxseen
for height,color in zip(heights,colors):
 ar=theory[height];ax.plot(ar[:,0],ar[:,1],color=color,ls='--',lw=1.7)
observed=max((row.get('tau') or 0 for row in records),default=0)
ax.set_xlim(0,observed*1.04 if observed else min(m.crossing(.005) for m in models.values())*.2)
handles=[Line2D([],[],color=c,marker=m,mfc='white',lw=1,label=r'$h_i/h_0=%.1f$'%h) for h,c,m in zip(heights,colors,markers)]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,1.0),ncol=3,columnspacing=.65,handlelength=1.2)
fig.text(.16,.78,r'BE, $\varphi=0.01$; all 25 seeds',va='bottom')
ax.axhline(.005,color='.5',ls=':',lw=.8);save(fig,'overview')
for height,color,marker in zip(heights,colors,markers):
 fig,ax=axes();ends=[]
 for cfg,style in zip([c for c in cases if c['initial']==height],['-','--','-.',':',(0,(5,1,1,1))]):ends.append(curve(ax,cfg,color,style,.7))
 ar=theory[height];ax.plot(ar[:,0],ar[:,1],color='k',ls='--',lw=2)
 ax.set_xlim(0,max(ends)*1.04 if max(ends)>0 else ar[-1,0]*.2)
 ax.axhline(.005,color='.5',ls=':',lw=.8)
 legend=[Line2D([],[],color=color,ls=ls,lw=1,label='Seed %d'%i) for i,ls in enumerate(['-','--','-.',':',(0,(5,1,1,1))],1)]
 legend.append(Line2D([],[],color='k',ls='--',lw=2,label='TCW + mode'))
 fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.55,1),ncol=3,columnspacing=.65,handlelength=1.2)
 fig.text(.16,.78,r'$h_i/h_0=%.1f$, $\varphi=0.01$'%(height),va='bottom');save(fig,'h%02d'%round(height*10))
panels=['overview']+['h%02d'%round(h*10) for h in heights]
sheet=Image.new('RGB',(3*960,2*720),'white')
for i,name in enumerate(panels):sheet.paste(Image.open(OUT/(name+'_preview.png')).convert('RGB'),((i%3)*960,(i//3)*720))
sheet.save(OUT/'board.png')
validation={name:item.get('validation.json',item.get('validation_failed.json',{})) for name,item in items.items() if item['route']['kind']=='validation'}
public=dict(updated_utc=snapshot['updated_utc'],parameters={k:snapshot['campaign'][k] for k in ['grid','noise_grid','cores_per_case','unit_time','stop']},records=records,validation=validation,curve_definition='Thin: all saved seeds; dashed: one registered TCW+mode per height, C_rms=3.5, white6, fixed 256² spectral sum. No ensemble mean or MFPT is claimed.')
(OUT/'data.json').write_text(json.dumps(public,indent=2),encoding='utf-8')
def fmt(v):return '—' if v is None or v=='' else ('%.5g'%v if isinstance(v,float) else str(v))
table=''.join('<tr>'+''.join('<td>'+html.escape(fmt(row.get(k)))+'</td>' for k in ['initial','replicate','server','job','node','status','tau','hmin','accepted','rejected','actual_wall_hours','planning_hours'])+'</tr>' for row in records)
passed=sum(v.get('passed',False) for v in validation.values());running=sum(r['status']=='RUNNING' for r in records);done=sum(r['status']=='THRESHOLD' for r in records)
panels_html=''.join('<article><h2>'+('全部厚度：全部25条轨迹' if name=='overview' else '初始厚度 '+str(int(name[1:])/10)+'：5条独立轨迹')+'</h2><a href="'+name+'.png"><img src="'+name+'_preview.png"></a><p><a href="'+name+'.pdf">PDF</a> · <a href="'+name+'.png">600 dpi PNG</a></p></article>' for name in panels)
page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>附壁液膜 BE 25轨迹看板</title><style>body{font-family:system-ui,"Microsoft YaHei",sans-serif;margin:0;background:#f4f6fa;color:#172436}header,main{max-width:1500px;margin:auto;padding:24px}.cards,.plots{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:16px}article,.card{background:white;border:1px solid #dce2eb;border-radius:12px;padding:18px}h1{font-size:28px}h2{font-size:17px}img{width:100%;height:auto}a{color:#185cac}.metric{font-size:30px;font-weight:650}.note{line-height:1.7;color:#526174}.table{overflow:auto;background:white;padding:16px;margin-top:20px}table{border-collapse:collapse;width:100%;font-size:13px;white-space:nowrap}td,th{padding:9px;text-align:left;border-bottom:1px solid #e2e7ee}button{padding:8px 15px;border-radius:6px;border:1px solid #abb6c7;background:white}@media(max-width:1000px){.plots{grid-template-columns:repeat(2,1fr)}}@media(max-width:650px){.cards,.plots{grid-template-columns:1fr}}@media(prefers-color-scheme:dark){body{background:#131c29;color:#edf1f6}article,.card,.table{background:#1c2939;border-color:#37465b}.note{color:#b5c0cf}a{color:#85b9ff}}</style><header><h1>附壁液膜 · 后向欧拉 BE · 25条独立轨迹</h1><p class="note">快照 UTC: STAMP。hᵢ/h₀=0.5–0.9，各5条；φ=0.01，无滑移，64²网格、white6噪声，4物理核/条。h₀为固定参考厚度。停止阈值hmin/h₀=0.005。</p><div class="cards"><div class="card">短测通过厚度<div class="metric">PASS / 5</div></div><div class="card">正在运行<div class="metric">RUN / 25</div></div><div class="card">实测到停止阈值<div class="metric">DONE / 25</div></div></div><p class="note">细线显示全部种子，不把删失轨迹计为破碎时间；虚线为各厚度唯一TCW+mode，Cᵣₘₛ=3.5。短测只筛查有限时间强误差，不代表首达统计收敛。末点空心圆=实测阈值，方形=尚在推进，叉=数值停止。绘图止于核验检查点。</p></header><main><div class="plots">PANELS</div><div class="table"><h2>逐条状态与资源预算</h2><p class="note">计划小时数按独立0.7试跑2小时54分钟乘TCW时间比例标定，用于资源安排；不含排队，不是MFPT或确定完工预报。实际吞吐、减步、节点负载与随机首达会改变预算。无保存检查点时显示空值。</p><p><a href="status.csv">状态CSV</a> · <a href="data.json">参数、验证与来源SHA</a></p><table><thead><tr><th>初厚</th><th>副本</th><th>服务器</th><th>作业</th><th>节点</th><th>状态</th><th>τ</th><th>hmin</th><th>接受</th><th>拒绝</th><th>实际小时</th><th>计划小时</th></tr></thead><tbody>TABLE</tbody></table></div></main></html>'''
for key,value in {'STAMP':snapshot['updated_utc'],'PASS':str(passed),'RUN':str(running),'DONE':str(done),'PANELS':panels_html,'TABLE':table}.items():page=page.replace(key,value)
(OUT/'index.html').write_text(page,encoding='utf-8')
print(json.dumps(dict(validation_passed=passed,running=running,threshold=done,snapshot=snapshot['updated_utc'],planned_hours_by_height={h:pilot_hours*models[h].crossing(.005)/models[.7].crossing(.005)*divisors.get(h,8)/8 for h in heights}),indent=2))
