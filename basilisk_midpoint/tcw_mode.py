"""Continuum 64-square white-noise TCW+mode; same normalization as midpoint.

No fitted seed dependence, no interpolation filter. C=5 is an envelope choice,
not an assertion that a single stochastic minimum equals the envelope.
"""
import numpy as np
from scipy.optimize import brentq

class TCWMode:
 def __init__(self,height,phi=.01,C=5,slip=0):
  self.height=height;self.C=C;n=64
  q=2*np.pi*np.fft.fftfreq(n)*n;q2=(q[:,None]**2+q[None,:]**2).ravel()[1:]
  unique,inverse=np.unique(q2,return_inverse=True)
  mobility=height**3+slip*height**2;A=13.159
  self.g=mobility*unique*(3*A/height**4-unique)
  self.Q=np.bincount(inverse,weights=2*phi*mobility*q2)
  self.omega=float(self.g.max());previous=1e-16
  for value in np.geomspace(1e-15,.1,500):
   if self.residual(previous)>0 and self.residual(value)<=0:
    self.switch=brentq(self.residual,previous,value,xtol=1e-16);break
   previous=value
  else:raise RuntimeError('TCW matching root not found')
  self.amplitude=C*np.sqrt(self.variance(self.switch)[0])
 def variance(self,t):
  rate=2*self.g;nonzero=abs(rate)>1e-15
  integral=np.full_like(rate,t);integral[nonzero]=np.expm1(rate[nonzero]*t)/rate[nonzero]
  modes=self.Q*integral
  return float(modes.sum()),float((2*self.g*modes+self.Q).sum())
 def residual(self,t):
  v,d=self.variance(t);return d/(2*v)-self.omega
 def predict(self,t):
  return self.height-self.C*np.sqrt(self.variance(t)[0]) if t<=self.switch else self.height-self.amplitude*np.exp(self.omega*(t-self.switch))
 def crossing(self,threshold):
  assert 0<=threshold<self.height
  if threshold>=self.predict(self.switch):return brentq(lambda t:self.predict(t)-threshold,0,self.switch,xtol=1e-16)
  return self.switch+np.log((self.height-threshold)/self.amplitude)/self.omega
