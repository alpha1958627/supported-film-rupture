"""Publish only a theory-vs-numerical mean rupture-time figure; inputs stay private."""
from pathlib import Path
import argparse,json,sys,os
import numpy as np

p=argparse.ArgumentParser();p.add_argument('--snapshot',required=True,help='Private verified Basilisk snapshot, not committed');args=p.parse_args()
repo=Path(__file__).resolve().parents[1];source=Path(args.snapshot);data=json.loads(source.read_text());rows=[r for r in data['records'] if r['kind']=='production'];assert len(rows)==25
sys.path.insert(0,str(repo/'basilisk_midpoint'));from tcw_mode import TCWMode
sys.path.insert(0,str(repo/'postprocess'));from academic_style import apply_academic_style
os.environ.setdefault('MPLCONFIGDIR',str(source.parent/'mpl_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

heights=[.5,.6,.7,.8,.9];means=[];audit=[]
for h in heights:
 group=[r for r in rows if r['latest']['config']['initial']==h];assert len(group)==5
 times=[]
 for r in group:
  d=r['latest'];e=d['event'];assert d['status']=='THRESHOLD' and d['config']['slip_ref']==0 and d['config']['phi_ref']==.01
  assert e['threshold']==.4 and e['hmin_before']>.4>=e['hmin_after'] and e['lower']<e['upper']==d['time']
  assert e['upper']-e['lower']<=h*1.01e-5/1e8
  times.append(e['upper'])
 means.append(float(np.mean(times)));audit.append(dict(initial=h,n=5,mean_tau=means[-1],events=times,theory_tau=TCWMode(h,C=5).crossing(.4)))
# Audit is saved next to the private snapshot, never in the published results folder.
(source.parent/'rupture_plot_private_audit.json').write_text(json.dumps(dict(snapshot_utc=data['utc'],groups=audit),indent=2))
apply_academic_style();fig=plt.figure(figsize=(8,6));ax=fig.add_axes([.16,.20,.80,.57]);x=np.linspace(.5,.9,161)
ax.plot(x,[TCWMode(h,C=5).crossing(.4) for h in x],color='k',lw=2,label=r'TCW + mode, $C_{\rm rms}=5$')
ax.plot(heights,means,color='#D62728',marker='o',lw=1.8,ms=6,label=r'Basilisk, mean ($n=5$)')
ax.set(xlabel=r'$h_i/h_0$',ylabel=r'Rupture time $\tau_r$',yscale='log',xlim=(.48,.92));ax.set_xticks(heights);ax.grid(False,which='both')
fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.55,.995),frameon=False)
fig.text(.55,.035,r'No slip, $\varphi=0.01$, $h_{\min}/h_0=0.4$',ha='center')
out=repo/'results';out.mkdir(exist_ok=True)
fig.savefig(out/'rupture_time_vs_h0.pdf');fig.savefig(out/'rupture_time_vs_h0.png',dpi=600);plt.close(fig)
print('Saved only PDF and PNG to results; five-seed event audit stays private')
