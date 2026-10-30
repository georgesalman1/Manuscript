import os, sys, json, numpy as np; ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,os.path.join(ROOT,'src'))
from mmg import *
out={}
P5=sine_profile(5.0)
print("== k=5 at ratio 0.02, T=30")
for sc in ['pa','cu','mms']:
    E=[run(n,sc,P5,ratio=0.02,T=30)['L2'] for n in [50,100,200,400]]
    print(sc,['%.2e'%e for e in E],['%.2f'%p for p in orders(E)[1:]]); out['k5_'+sc]=E
M=multiscale_profile(); N2=[50,100,200,400,800,1600]
print("== multiscale PA T=5 vs T=10, sum u conservation")
for T in [5,10]:
    R=[run(n,'pa',M,bc='periodic',T=T,ratio=0.05) for n in N2]
    E=[r['L2'] for r in R]; print(T,['%.2e'%e for e in E],['%.2f'%p for p in orders(E)[1:]],'drift max %.1e'%max(r['drift'] for r in R))
    out[f'ms_pa_T{T}']=E
    su=[abs(r['u'].sum()-r['ue'].sum())/abs(r['ue'].sum()) for r in R]; print('  sum-u rel change',['%.1e'%x for x in su])
P3=sine_profile(3.0)
R=[run(n,'pa',P3,bc='periodic',L=2*np.pi,T=40,ratio=0.04) for n in [50,100,200,400]]
print('per k3 pa sum-u rel change',['%.1e'%(abs(r['u'].sum()-r['ue'].sum())/abs(r['ue'].sum())) for r in R], 'drift',['%.1e'%r['drift'] for r in R])
json.dump(out,open(os.path.join(ROOT,'results','exp7.json'),'w'))
