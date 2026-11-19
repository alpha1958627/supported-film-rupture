"""C=5 review using the exact supplied-figure Basilisk rows; preserve prior versions."""
from pathlib import Path
import argparse,json,csv,io,hashlib,subprocess,base64,sys,os
import numpy as np

p=argparse.ArgumentParser();p.add_argument('--workspace',required=True);p.add_argument('--source-root',default='F:/Predict_rupture/spinodal_comparison_20261110');args=p.parse_args()
ROOT=Path(args.workspace).resolve();REPO=Path(__file__).resolve().parent;SOURCE=Path(args.source_root)
OUT=REPO/'dashboard/review_C5_image_reference';OUT.mkdir(exist_ok=True)
files=dict(figure_csv=SOURCE/'tcw_extremes/event_times_2024PRF.csv',figure_script=SOURCE/'tcw_extremes/analyze_tcw_prf.py',original_preview=SOURCE/'tcw_extremes/event_times_2024PRF_preview.png',individual=SOURCE/'fit_h05/individual_events.csv',crossing_script=SOURCE/'fit_h05/plot_fit_h05.py',h05_cases=SOURCE/'fit_h05/basilisk_h05_cases.csv')
def csvread(path):return list(csv.DictReader(path.open(encoding='utf-8-sig')))
def csvwrite(path,rows):
 stream=io.StringIO();writer=csv.DictWriter(stream,fieldnames=list(dict.fromkeys(k for row in rows for k in row)),lineterminator='\n');writer.writeheader();writer.writerows(rows);path.write_bytes(stream.getvalue().encode('utf-8-sig'))
original=csvread(files['figure_csv']);bm=[v for v in original if v['series']=='Basilisk'];individual=[v for v in csvread(files['individual']) if v['series']=='Basilisk'];image_model=[v for v in original if v['series']=='TCW 5']
assert len(bm)==5 and len(individual)==25 and sorted(float(v['h0']) for v in bm)==[.5,.6,.7,.8,.9]
for row in bm:
 values=[float(v['tau']) for v in individual if float(v['h0'])==float(row['h0'])]
 assert len(values)==5 and int(float(row['n']))==5 and np.isclose(np.mean(values),float(row['tau']),rtol=1e-10,atol=1e-15)
 assert row['criterion']=='first_recorded_Hmin_0.4'
csvwrite(OUT/'Basilisk_exact_figure_means.csv',bm);csvwrite(OUT/'Basilisk_original_interpolated_events.csv',individual)
csvwrite(OUT/'original_figure_TCW_C5.csv',image_model)
snapshot=json.loads((ROOT/'snapshot.json').read_text());items=[v for g in snapshot['groups'] for v in g.get('items',[]) if v['route']['kind']=='production'];assert len(items)==25
records=[]
for item in items:
 d=item['latest'];records.append(dict(name=item['route']['name'],remote=item['route']['remote'],height=d['config']['initial'],seed=d['config']['seed'],accepted=d['accepted'],time=d['time'],history_prefix_sha256=item['history_prefix_sha256']))
remote=r'''
from pathlib import Path
import json,hashlib,io
import numpy as np
items=ITEMS;rows=[]
for v in items:
 with (Path(v['remote'])/'run/history.csv').open('rb') as f:raw=b''.join(f.readline() for _ in range(v['accepted']+1))
 assert hashlib.sha256(raw).hexdigest()==v['history_prefix_sha256']
 data=np.loadtxt(io.BytesIO(raw),delimiter=',',skiprows=1,usecols=(0,2),ndmin=2)
 assert len(data)==v['accepted'] and data[-1,0]==v['time']
 threshold=.4*v['height'];hit=np.flatnonzero(data[:,1]<=threshold)
 row=dict(name=v['name'],height=v['height'],seed=v['seed'],threshold=threshold,reached=bool(len(hit)),observed_tau=v['time'],history_prefix_sha256=v['history_prefix_sha256'])
 if len(hit):
  i=int(hit[0]);before=data[i-1] if i else [0,v['height']]
  row.update(tau_lower=float(before[0]),tau_upper=float(data[i,0]),h_before=float(before[1]),h_after=float(data[i,1]),accepted_step=i+1)
 rows.append(row)
print(json.dumps(rows))
'''
cachepath=ROOT/'BE_first_local_H04_review.json';cache=json.loads(cachepath.read_text()) if cachepath.exists() else []
lookup={v['name']:v for v in cache}
valid=len(cache)==25 and all(lookup.get(v['name'],{}).get('history_prefix_sha256')==v['history_prefix_sha256'] and lookup[v['name']]['observed_tau']==v['time'] for v in records)
if not valid:
 code=remote.replace('ITEMS',repr(records));command='echo '+base64.b64encode(code.encode()).decode()+' | base64 -d | '+items[0]['route']['python']+' -'
 r=subprocess.run(['ssh','-o','BatchMode=yes',items[0]['route']['host'],command],capture_output=True,text=True,timeout=240)
 if r.returncode:raise RuntimeError(r.stderr)
 cache=json.loads(r.stdout);cachepath.write_text(json.dumps(cache,indent=2),encoding='utf-8')
csvwrite(OUT/'BE_first_local_H04.csv',cache)
absolute=json.loads((ROOT/'crossings_h04_raw.json').read_text())['BE']
assert len(absolute)==25 and all(next(v for v in absolute if v['name']==r['name'])['history_prefix_sha256']==r['history_prefix_sha256'] for r in records)
sys.path.insert(0,str(REPO/'be'));from tcw_mode import TCWMode
heights=[.5,.6,.7,.8,.9];models={h:TCWMode(h,C=5) for h in heights};summary=[]
for h,m in models.items():
 reached=[v['tau_upper'] for v in cache if v['height']==h and v['reached']];fixed=[v['tau_upper'] for v in absolute if v['height']==h and v['reached']]
 target=.4*h;t=m.crossing(target);assert abs(m.predict(t)-target)<1e-8
 summary.append(dict(height=h,common_threshold_h_over_reference=target,Basilisk_image_tau=float(next(v for v in bm if float(v['h0'])==h)['tau']),BE_local_H04_reached=len(reached),BE_local_H04_mean=float(np.mean(reached)) if len(reached)==5 else None,TCW_C5_current_parameters_local_H04=t,BE_absolute_04_reached=len(fixed),BE_absolute_04_mean=float(np.mean(fixed)) if len(fixed)==5 else None,TCW_C5_current_parameters_absolute_04=m.crossing(.4)))
provenance=dict(snapshot_utc=snapshot['updated_utc'],C_rms=5,phi=.01,original_figure_definition='H=h/h_initial; Hmin=.4. Original native initial H=1. Common u=h/h_reference therefore has threshold u=.4*height.',source_files={k:dict(path=str(v),sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k,v in files.items()},Basilisk_rows_unchanged=bm,original_crossing_method='First saved native H<=.4, linear interpolation of the preceding sampled segment; mean of five paths. These are supplied-figure estimates, not exactly resolved stochastic first-passage times.',theory_used='Current no-slip m=h^3, phi=.01, white6 sinc^4 transfer, fixed 256² spectrum; C changed to5. The supplied figure original theory instead used slip=.2 and a 64² continuum spectrum without this interpolation transfer.',original_h05_source='Nov10 supplementary five paths, including seed4 retry; .6-.9 from Oct25 formal batch.',summary=summary)
(OUT/'provenance_and_summary.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
os.environ['MPLCONFIGDIR']=str(ROOT/'mpl_cache');sys.path.insert(0,str(REPO/'postprocess'))
from academic_style import apply_academic_style
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
apply_academic_style()
def newfig():
 fig=plt.figure(figsize=(8,6));ax=fig.add_axes([.16,.19,.81,.52]);ax.set(xlabel=r'$h_i/h_0$',ylabel=r'Dimensionless time $\tau$',yscale='log',xlim=(.46,.94));ax.set_xticks(heights);ax.grid(False,which='both');return fig,ax
def save(fig,name):
 fig.savefig(OUT/(name+'.pdf'));fig.savefig(OUT/(name+'.png'),dpi=600);fig.savefig(OUT/(name+'_preview.png'),dpi=120);plt.close(fig)
def observations(ax,rows,color):
 for h in heights:
  hits=[v for v in rows if v['height']==h and v['reached']];incomplete=[v for v in rows if v['height']==h and not v['reached']]
  for i,v in enumerate(hits):ax.plot(h+(i-2)*.003,v['tau_upper'],'o',mfc='white',mec=color,ms=5)
  if len(hits)==5:
   values=np.array([v['tau_upper'] for v in hits]);ax.errorbar(h,values.mean(),yerr=values.std(ddof=1),color=color,fmt='s',ms=5,capsize=4,lw=1.3)
  for i,v in enumerate(incomplete):
   x=h+(i-2)*.003;y=v['observed_tau'];ax.plot(x,y,'^',mfc='white',mec=color,ms=5);ax.annotate('',xy=(x,y*1.4),xytext=(x,y),arrowprops=dict(arrowstyle='->',color=color,lw=.8))
grid=np.linspace(.5,.9,25)
fig,ax=newfig();ax.plot(grid,[TCWMode(h,C=5).crossing(.4*h) for h in grid],color='k',lw=2)
ax.plot([float(v['h0']) for v in bm],[float(v['tau']) for v in bm],color='#D62728',marker='o',lw=1.5,ms=6)
observations(ax,cache,'#1F77B4')
handles=[Line2D([],[],color='k',lw=2,label=r'TCW + mode, $C_{\rm rms}=5$'),Line2D([],[],color='#D62728',marker='o',lw=1.5,label='Basilisk: supplied figure'),Line2D([],[],color='#1F77B4',marker='o',mfc='white',ls='none',label='BE: first crossing'),Line2D([],[],color='#1F77B4',marker='^',mfc='white',ls='none',label='BE: incomplete')]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.995),labelspacing=.15);fig.text(.55,.035,r'Common criterion: $h_{\min}/h_i=0.4$.',ha='center');save(fig,'C5_common_initial_fraction_04')
fig,ax=newfig();ax.plot(grid,[TCWMode(h,C=5).crossing(.4) for h in grid],color='k',lw=2);observations(ax,absolute,'#1F77B4')
fig.legend(handles=[handles[0],handles[2],handles[3]],loc='upper center',bbox_to_anchor=(.55,.995),labelspacing=.15);fig.text(.55,.035,r'Previous criterion: $h_{\min}/h_0=0.4$.',ha='center');save(fig,'C5_previous_absolute_04')
fig,ax=newfig();ax.plot([float(v['h0']) for v in image_model],[float(v['tau']) for v in image_model],color='k',ls='-.',lw=2);ax.plot([float(v['h0']) for v in bm],[float(v['tau']) for v in bm],color='#D62728',marker='o',lw=1.5,ms=6)
fig.legend(handles=[Line2D([],[],color='k',ls='-.',lw=2,label=r'Image TCW + mode, $C_{\rm rms}=5$'),handles[1]],loc='upper center',bbox_to_anchor=(.55,.995));fig.text(.55,.035,r'Original figure: $H_{\min}=0.4$.',ha='center');save(fig,'source_figure_C5_check')
lines=['# C_rms=5与所附原图Basilisk数据核对版','', 'BE沿用快照UTC '+snapshot['updated_utc']+'。原C=3.5及绝对阈值0.3/0.4版本保持不变。','', '## 原图红点，直接读取CSV而非图片取点','', '|初厚|原图Basilisk平均tau|路径数|对应当前厚度阈值|','|---|---:|---:|---:|']
for row in bm:lines.append('|%.1f|%.12g|5|%.2f|'%(float(row['h0']),float(row['tau']),.4*float(row['h0'])))
lines+=['', '原图H以各自初始厚度归一：H=h/h_i，所以H=.4对应当前h/h0=.4(h_i/h0)，并非固定h/h0=.4。原图Basilisk时间用相邻稀疏观测线性插值，再对5条取平均；本版按原样恢复，不将其重称为精确首达或严格括区。0.5来自11月10日补充批次，0.6–0.9来自10月25日正式批次。', '', '## 同判据：hmin/hi=.4，C=5','', '|初厚|BE到达数|BE平均tau（全部5条到达才统计）|当前无滑移TCW+mode C5|','|---|---:|---:|---:|']
f=lambda x:'—' if x is None else '%.9g'%x
for v in summary:lines.append('|%.1f|%d/5|%s|%s|'%(v['height'],v['BE_local_H04_reached'],f(v['BE_local_H04_mean']),f(v['TCW_C5_current_parameters_local_H04'])))
lines+=['', '主图黑线保持当前BE的无滑移、white6传递和固定256²理论求和，只把C改为5。所附图原理论有滑移.2、64²连续Fourier谱，故source_figure_C5_check另按原图CSV原样重画，不能把两套黑线视为仅C不同。', '', '本目录为独立核对版，不覆盖原0.3/0.4参考，也不修改模拟运行参数。BE蓝点是每条首次接受步下穿，方块及误差棒为5条均值±样本SD，上箭头表示未到达的观察下界。', '', '[共同判据图](C5_common_initial_fraction_04.png) · [原固定0.4判据，仅BE与C5理论](C5_previous_absolute_04.png) · [所附图C5与红点原样核对](source_figure_C5_check.png) · [来源与统计](provenance_and_summary.json)']
(OUT/'README.md').write_text('\n'.join(lines),encoding='utf-8')
page=REPO/'dashboard/index.html';html=page.read_text(encoding='utf-8')
section='<section id="C5-image-review"><h2>C_rms=5：所附原图Basilisk数据核对版</h2><p class="note">共同判据hmin/hi=.4，对应当前参考厚度阈值.4×hi/h0。Basilisk红点直接恢复原CSV的5条均值与插值定义；黑线保持当前无滑移BE参数，只将C改5。</p><a href="review_C5_image_reference/C5_common_initial_fraction_04.png"><img style="max-width:960px" src="review_C5_image_reference/C5_common_initial_fraction_04_preview.png"></a><p><a href="review_C5_image_reference/README.md">数据和定义核对</a> · <a href="review_C5_image_reference/source_figure_C5_check.png">原图C5曲线与红点原样核对</a> · <a href="review_C5_image_reference/C5_previous_absolute_04.png">原固定绝对0.4判据C5版</a></p></section>'
if 'id="C5-image-review"' in html:
 start=html.index('<section id="C5-image-review">');end=html.index('</section>',start)+len('</section>');html=html[:start]+html[end:]
page.write_text(html.replace('</main>',section+'</main>'),encoding='utf-8')
print(json.dumps(summary,indent=2))
