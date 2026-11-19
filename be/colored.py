"""JCP discrete SE covariance with replayable stationary OU endpoints."""
import math
import numpy as np

class OUTape:
    def __init__(self,seed,n,base_dt=.003125,depth=20,tau=.01,sigma=.05,length=128.):
        self.seed=seed;self.n=n;self.base_dt=base_dt;self.depth=depth;self.tau=tau;self.sigma=sigma
        self.width=1<<depth;self.quantum=base_dt/self.width;self.cursor=0;self.cell=0;self.cache={}
        q=2*np.pi*np.fft.fftfreq(n)*n/length;w=np.exp(-.125*sigma*sigma*q*q)
        self.filter=w[:,None]*w;self.filter[0,0]=0.
        self.norm=math.sqrt(float(np.sum(w*w))**2-1)/n
        self.left=self.normal(0,0,0,91)
    def ticks(self,t):return int(round(t/self.quantum))
    def normal(self,cell,level,index,kind):
        rng=np.random.default_rng(np.random.SeedSequence([self.seed,cell&0xffffffff,cell>>32,level,index,kind]))
        a=rng.standard_normal((2,self.n,self.n))
        return np.fft.ifft2(np.fft.fft2(a,axes=(-2,-1))*self.filter,axes=(-2,-1)).real/self.norm
    def node(self,level,index):
        key=(level,index)
        if key in self.cache:return self.cache[key]
        if level==0:
            a=math.exp(-self.base_dt/self.tau)
            self.cache[key]=(self.left,a*self.left+math.sqrt(-math.expm1(-2*self.base_dt/self.tau))*self.normal(self.cell,0,0,92))
        else:
            left,right=self.node(level-1,index//2)
            leftkey=(level,2*(index//2));rightkey=(level,2*(index//2)+1)
            if leftkey not in self.cache and rightkey not in self.cache:
                duration=self.base_dt/(1<<(level-1));b=math.exp(-duration/(2*self.tau))
                mid=b/(1+b*b)*(left+right)+math.sqrt(-math.expm1(-duration/self.tau)/(1+b*b))*self.normal(self.cell,level-1,index//2,93)
                self.cache[leftkey]=(left,mid);self.cache[rightkey]=(mid,right)
        return self.cache[key]
    def value(self,start,end):
        assert start==self.cursor and start//self.width==self.cell
        size=end-start;assert size>0 and size&(size-1)==0 and start%size==0
        level=self.depth-(size.bit_length()-1);index=(start%self.width)//size
        return self.node(level,index)[1]
    def commit(self,end):
        assert end>=self.cursor;self.cursor=end
        if end==(self.cell+1)*self.width:
            self.left=self.node(0,0)[1];self.cell+=1;self.cache={}
        else:
            local=end-self.cell*self.width
            self.cache={key:value for key,value in self.cache.items() if (key[1]+1)*(self.width>>key[0])>local}

class FaceInterpolator:
    def __init__(self,n,coarse):
        self.n=n;self.coarse=coarse;self.maps=[]
        for axis in (0,1):
            gx=np.arange(n)*coarse/n if axis==0 else (np.arange(n)+.5)*coarse/n-.5
            gy=(np.arange(n)+.5)*coarse/n-.5 if axis==0 else np.arange(n)*coarse/n
            ix=np.floor(gx).astype(int);iy=np.floor(gy).astype(int)
            self.maps.append((ix%coarse,iy%coarse,(gx-ix)[:,None],(gy-iy)[None,:]))
    def __call__(self,a):
        assert a.shape==(2,self.coarse,self.coarse)
        out=np.empty((2,self.n,self.n))
        for axis,(ix,iy,fx,fy) in enumerate(self.maps):
            b=a[axis];jx=(ix+1)%self.coarse;jy=(iy+1)%self.coarse
            out[axis]=(1-fx)*(1-fy)*b[ix[:,None],iy]+fx*(1-fy)*b[jx[:,None],iy]+(1-fx)*fy*b[ix[:,None],jy]+fx*fy*b[jx[:,None],jy]
        return np.ascontiguousarray(out)
