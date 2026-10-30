import os, sys, numpy as np; ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,os.path.join(ROOT,'src'))
from mmg import *
def crit(sc,P,n=100,hi=0.5,T=40,**kw):
    lo=0.0
    for _ in range(12):
        mid=(lo+hi)/2
        if run(n,sc,P,ratio=mid,T=T,**kw)['blowup']: hi=mid
        else: lo=mid
    return lo
print("== stability limit vs k (n=100)")
for k in [1.05,1.25,1.5,2,3,5,10]:
    P=sine_profile(k); o=run(100,'pa',P,T=0.01)
    cp=crit('pa',P); cc=crit('cu',P)
    s=k+np.sin(np.arange(1,101)/100); cmax=(1/s+np.sqrt(s)).max()
    print(k,'Fr %.2f'%o['Fr_max'],'pa %.4f (CFL %.3f)'%(cp,cp*cmax),'cu %.4f (CFL %.3f)'%(cc,cc*cmax))
