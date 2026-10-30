import os, sys, json, numpy as np; ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,os.path.join(ROOT,'src'))
from mmg import *
P=sine_profile(2.0); N=[50,100,200,400]; out={}
print("== mutation (fault-injection) tests, T=2 ratio 0.02")
for mut in [None,'grav2','bedsign','adv2dx','hcentral']:
    for sc in ['pa','mms']:
        res=[run(n,sc,P,mutation=mut) for n in N]
        if any(r['blowup'] for r in res): print(mut,sc,'BLOWUP'); out[f'mut_{mut}_{sc}']='blowup'; continue
        E=[r['L2'] for r in res]; p=orders(E)[1:]
        print(mut,sc,['%.2e'%e for e in E],['%.2f'%x for x in p]); out[f'mut_{mut}_{sc}']=(E,p)
print("== zero-crossing MMS fields")
cases={'converging u=0.5cos(pi x)':(lambda x:0.5*np.cos(np.pi*x), lambda x:-0.5*np.pi*np.sin(np.pi*x)),
       'diverging u=-0.5cos(pi x)':(lambda x:-0.5*np.cos(np.pi*x), lambda x:0.5*np.pi*np.sin(np.pi*x))}
for name,(f,df) in cases.items():
    for d in ['fwd','adaptive']:
        res=[run(n,'mms',P,direction=d,uexact=f,duexact=df,T=20,ratio=0.5) for n in N]
        if any(r['blowup'] for r in res): print(name,d,'BLOWUP at t=%.2f'%res[0]['t']); out[f'zc_{name}_{d}']='blowup'; continue
        E=[r['L2'] for r in res]; p=orders(E)[1:]
        print(name,d,['%.2e'%e for e in E],['%.2f'%x for x in p]); out[f'zc_{name}_{d}']=(E,p)
json.dump(out,open(os.path.join(ROOT,'results','exp4.json'),'w'),default=float)
