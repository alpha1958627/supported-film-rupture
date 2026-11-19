import math
import numpy as np

class BrownianTape:
    def __init__(self,seed,shape,base_dt=.001,depth=12):
        assert 0<=int(seed)<2**64 and 1<=depth<=22
        self.seed=int(seed);self.shape=tuple(shape);self.base_dt=float(base_dt);self.depth=int(depth)
        self.width=1<<self.depth
        self.quantum=self.base_dt/self.width
        self.cursor=0
        self.cache={}

    def ticks(self,time):
        return int(round(float(time)/self.quantum))

    def _normal(self,cell,level,node,kind):
        assert 0<=cell<2**128
        entropy=[self.seed&0xffffffff,self.seed>>32]+[(cell>>(32*i))&0xffffffff for i in range(4)]+[level,node,kind]
        return np.random.default_rng(np.random.SeedSequence(entropy)).standard_normal(self.shape)

    def node(self,cell,level,index):
        key=(cell,level,index)
        if key in self.cache:
            return self.cache[key]
        if level==0:
            self.cache[key]=math.sqrt(self.base_dt)*self._normal(cell,0,0,1)
        else:
            parent=self.node(cell,level-1,index//2)
            left=(cell,level,2*(index//2));right=(cell,level,2*(index//2)+1)
            if left in self.cache:
                self.cache[right]=parent-self.cache[left]
            elif right in self.cache:
                self.cache[left]=parent-self.cache[right]
            else:
                duration=self.base_dt/(1<<(level-1))
                self.cache[left]=.5*parent+math.sqrt(duration/4)*self._normal(cell,level-1,index//2,2)
                self.cache[right]=parent-self.cache[left]
        return self.cache[key]

    def _inside(self,cell,a,b,level=0,index=0,start=0,width=None):
        width=self.width if width is None else width
        if a==start and b==start+width:
            return self.node(cell,level,index)
        half=width//2
        assert half>0
        middle=start+half
        if b<=middle:
            return self._inside(cell,a,b,level+1,2*index,start,half)
        if a>=middle:
            return self._inside(cell,a,b,level+1,2*index+1,middle,half)
        return (self._inside(cell,a,middle,level+1,2*index,start,half)+
                self._inside(cell,middle,b,level+1,2*index+1,middle,half))

    def increment(self,start,end):
        assert self.cursor<=start<end
        result=np.zeros(self.shape)
        for cell in range(start//self.width,(end-1)//self.width+1):
            a=max(start-cell*self.width,0);b=min(end-cell*self.width,self.width)
            result+=self._inside(cell,a,b)
        return result

    def commit(self,end):
        assert end>=self.cursor
        self.cursor=int(end)
        self.cache={key:value for key,value in self.cache.items() if (key[0]+1)*self.width>self.cursor}

