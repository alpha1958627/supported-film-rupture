"""Adaptive no-repulsion first-passage test. Rejected trials reuse the noise path."""
from pathlib import Path
import os,sys,json,time,csv,ctypes,pickle,hashlib,signal,math
import numpy as np
from noise import BrownianTape
from colored import OUTape,FaceInterpolator

root=Path(__file__).resolve().parent
_BOUND_CPUS=None
def atomic_json(path,obj):
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(obj,indent=2));os.replace(temp,path)
def affinity():
    global _BOUND_CPUS
    if _BOUND_CPUS is not None:return _BOUND_CPUS
    cores={}
    for c in sorted(os.sched_getaffinity(0)):
        p=Path(f'/sys/devices/system/cpu/cpu{c}/topology');cores.setdefault(((p/'physical_package_id').read_text(),(p/'core_id').read_text()),c)
    cpus=list(cores.values())[:4];assert len(cpus)==4
    os.environ.update(OMP_NUM_THREADS='4',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='TRUE',OMP_PLACES=','.join('{'+str(c)+'}' for c in cpus))
    _BOUND_CPUS=cpus
    return cpus

class Simulation:
    def __init__(self,cfg,out):
        self.cfg=cfg;self.out=out;out.mkdir(parents=True,exist_ok=True);self.n=1<<cfg['grid_level']
        self.cpus=affinity();self.lib=ctypes.CDLL(str(root/'wall_native.so'));ptr=ctypes.POINTER(ctypes.c_double)
        self.lib.wall_init.argtypes=[ctypes.c_int,ctypes.c_int];self.lib.wall_step.argtypes=[ptr,ptr,ctypes.c_double,ptr,ptr]
        self.ptr=lambda a:a.ctypes.data_as(ptr)
        assert self.lib.wall_init(cfg['grid_level'],2 if cfg['noise']=='ouse' else 0)==0
        assert cfg['noise']=='white'
        self.lib.wall_set_phi.argtypes=[ctypes.c_double];self.lib.wall_set_phi(cfg['phi'])
        self.cap=2.1e-9
        coarse=self.n if cfg['noise']=='ouse' else 1<<cfg['noise_level']
        self.noise=OUTape(cfg['seed'],coarse) if cfg['noise']=='ouse' else BrownianTape(cfg['seed'],(2,coarse,coarse),base_dt=self.cap,depth=22)
        self.interp=FaceInterpolator(self.n,coarse);self.coarse_dx=1/coarse
        self.h=np.full((self.n,self.n),cfg['initial']);self.cursor=0;self.max_proposal=self.noise.width//cfg['dt_divisor'];self.proposal=self.max_proposal;self.steps=0;self.rejected=0;self.streak=0
        self.status='RUNNING';self.event=None;self.wall_total=0.;self.actual_wall_total=0.;self.previous=None
        self.sources={n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in ('wall_native.c','wall_native.so','worker.py','noise.py','colored.py')}
        self.history=out/'history.csv';self.initial_mean=cfg['initial']
        if (out/'checkpoint.npz').exists():self.restore()
        else:
            with self.history.open('w',newline='') as f:csv.writer(f).writerow(['time','dt','hmin','hmax','relative_mass','full_residual','outer_iterations','mg_cycles','accepted','rejected','wall_seconds'])
        atomic_json(out/'allocation.json',dict(job=os.getenv('SLURM_JOB_ID'),host=os.uname().nodename,physical_cpus=self.cpus))
    @property
    def time(self):return self.cursor*self.noise.quantum
    def save(self):
        state={k:v for k,v in self.__dict__.items() if k not in ('lib','ptr','interp','h','out','history','cpus')}
        state['history_sha256']=hashlib.sha256(self.history.read_bytes()).hexdigest()
        payload=pickle.dumps(state,protocol=4);digest=hashlib.sha256(self.h.tobytes()+payload).hexdigest()
        temp=self.out/'checkpoint.tmp.npz'
        np.savez_compressed(temp,h=self.h,payload=np.frombuffer(payload,dtype=np.uint8),sha256=digest)
        os.replace(temp,self.out/'checkpoint.npz')
        atomic_json(self.out/'latest.json',dict(status=self.status,time=self.time,hmin=float(self.h.min()),hmax=float(self.h.max()),relative_mass=float(self.h.mean()/self.initial_mean-1),accepted=self.steps,rejected=self.rejected,event=self.event,wall_seconds=self.wall_total,actual_wall_seconds=self.actual_wall_total,updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),config=self.cfg,source_sha256=self.sources['wall_native.c'],checkpoint_sha256=hashlib.sha256((self.out/'checkpoint.npz').read_bytes()).hexdigest()))
    def restore(self):
        with np.load(self.out/'checkpoint.npz',allow_pickle=False) as z:
            h=z['h'].copy();payload=z['payload'].tobytes();assert hashlib.sha256(h.tobytes()+payload).hexdigest()==str(z['sha256'])
        state=pickle.loads(payload);assert state['cfg']==self.cfg and state['sources']==self.sources
        self.__dict__.update(state);self.h=h;self.out=Path(self.out);self.history=self.out/'history.csv'
        # Ignore at most an uncheckpointed tail, retaining it for audit.
        lines=self.history.read_bytes().splitlines(keepends=True);keep=self.steps+1
        if len(lines)>keep:
            (self.out/('uncheckpointed_history_'+str(int(time.time()*1e9))+'.csv')).write_bytes(b''.join(lines))
            self.history.write_bytes(b''.join(lines[:keep]))
        assert len(self.history.read_bytes().splitlines())==keep
        assert hashlib.sha256(self.history.read_bytes()).hexdigest()==state['history_sha256']
    def attempt(self,ticks):
        dt=ticks*self.noise.quantum
        if self.cfg['noise']=='ouse':raw=self.noise.value(self.cursor,self.cursor+ticks)
        else:raw=self.noise.increment(self.cursor,self.cursor+ticks)/(self.coarse_dx*dt)
        drive=self.interp(raw);h=np.empty_like(self.h);d=np.zeros(8)
        code=self.lib.wall_step(self.ptr(self.h),self.ptr(drive),dt,self.ptr(h),self.ptr(d))
        return code,h,d
    def advance(self):
        ticks=min(self.proposal,self.noise.width-(self.cursor%self.noise.width))
        while self.cursor%ticks:ticks//=2
        event_ticks=max(1,self.noise.ticks(1e-6/805306368.));why=''
        while True:
            started=time.monotonic();code,h,d=self.attempt(ticks);self.wall_total+=time.monotonic()-started
            valid=code==0 and np.isfinite(h).all() and np.isfinite(d).all() and h.min()>0 and d[3]<1e-8 and abs(h.mean()/self.initial_mean-1)<1e-9
            crossing=valid and h.min()<=self.cfg['threshold']
            if valid and not(crossing and ticks>event_ticks):break
            self.rejected+=1;self.streak=0
            why='event refinement' if crossing else 'positivity/nonlinear/residual/mass rejection'
            if ticks==1:self.status='NUMERICAL_FAILURE';raise RuntimeError(why)
            ticks//=2
        before=self.time;oldmin=float(self.h.min());self.previous=self.h.copy() if crossing else None
        self.h=h;self.cursor+=ticks;self.noise.commit(self.cursor);self.steps+=1
        dt=ticks*self.noise.quantum
        if crossing:
            self.status='THRESHOLD';self.event=dict(lower=before,upper=self.time,hmin_before=oldmin,hmin_after=float(h.min()),threshold=self.cfg['threshold'])
        self.proposal=ticks;self.streak+=1
        if self.streak>=20:self.proposal=min(self.max_proposal,2*ticks);self.streak=0
        row=[self.time,dt,float(h.min()),float(h.max()),float(h.mean()/self.initial_mean-1),float(d[3]),int(d[6]),int(d[7]),self.steps,self.rejected,self.wall_total]
        with self.history.open('a',newline='') as f:csv.writer(f).writerow(row)
        return row

def main():
    cfg=json.loads((Path.cwd()/'config.json').read_text());sim=Simulation(cfg,Path.cwd()/'run')
    stop=[False]
    def signal_stop(*args):stop[0]=True
    signal.signal(signal.SIGUSR1,signal_stop);signal.signal(signal.SIGTERM,signal_stop)
    started=lastsave=time.monotonic();lastprint=started;last_iteration=started
    try:
        while sim.status=='RUNNING' and not stop[0]:
            row=sim.advance();now=time.monotonic();sim.actual_wall_total+=now-last_iteration;last_iteration=now
            if now-lastsave>=300:sim.save();lastsave=now
            if now-lastprint>=60:print(json.dumps(dict(time=sim.time,hmin=row[2],dt=row[1],accepted=sim.steps,rejected=sim.rejected,actual_wall_seconds=sim.actual_wall_total,tau_per_actual_hour=sim.time*3600/max(sim.actual_wall_total,1e-9))),flush=True);lastprint=now
            if now-started>71.5*3600:stop[0]=True
        sim.save()
        if sim.status=='THRESHOLD':return 0
        if stop[0]:return 75
        return 82
    except Exception as e:
        sim.status='NUMERICAL_FAILURE';sim.save();atomic_json(sim.out/'FAILED.json',dict(error=repr(e),time=sim.time));raise

if __name__=='__main__':sys.exit(main())
