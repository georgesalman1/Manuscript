import os, sys, json, numpy as np; ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,os.path.join(ROOT,'src'))
from mmg import *
N=[50,100,200,400]; out={}
def conv(sc,P,Ns=N,**kw):
    res=[run(n,sc,P,**kw) for n in Ns]
    if any(r['blowup'] for r in res): return None,res
    E=[r['L2'] for r in res]; return E,res
print("== k sweep (steady state)")
for k in [1.05,1.1,1.25,1.5,2,3,5]:
    P=sine_profile(k); row={}
    for sc,ratio in [('pa',0.05 if k<=3 else 0.02),('cu',0.05 if k<=3 else 0.02),('mms',0.5)]:
        E,res=conv(sc,P,T=30,ratio=ratio)
        row[sc]=(E,orders(E)[1:] if E else None)
    Fr=res[0]['Fr_max']
    print(k,'Fr %.2f'%Fr,{s:(['%.2e'%e for e in v[0]],['%.2f'%p for p in v[1]]) if v[0] else 'X' for s,v in row.items()})
    out[f'k{k}']=dict(Fr=Fr,**{s:v for s,v in row.items()})
print("== rightward flow u>0 (C=+1)")
P=sine_profile(2.0)
for sc,ratio in [('pa',0.05),('cu',0.05),('mms',0.5)]:
    for d in ['fwd','bwd']:
        E,res=conv(sc,P,C=+1.0,direction=d,T=20,ratio=ratio)
        print(sc,d, (['%.2e'%e for e in E],['%.2f'%p for p in orders(E)[1:]]) if E else 'BLOWUP t=%.2f'%res[0].get('t',-1))
        out[f'right_{sc}_{d}']=(E,orders(E)[1:]) if E else 'blowup'
print("== periodic, k+sin(x), k=3 on (0,2pi], T=40")
P3=sine_profile(3.0)
for sc,ratio in [('pa',0.04),('cu',0.04)]:
    E,res=conv(sc,P3,bc='periodic',L=2*np.pi,T=40,ratio=ratio)
    print(' mean flux',['%.4f'%np.mean(r['u']*r['h']) for r in res])
    print(sc,(['%.2e'%e for e in E],['%.2f'%p for p in orders(E)[1:]]) if E else 'X', 'drift',['%.1e'%r.get('drift',float('nan')) for r in res])
    out[f'per_{sc}']=(E,orders(E)[1:] if E else None,[r['drift'] for r in res])
print("== multiscale periodic (0,1]")
M=multiscale_profile(); N2=[50,100,200,400,800,1600]
for sc,ratio in [('pa',0.05),('cu',0.05)]:
    E,res=conv(sc,M,Ns=N2,bc='periodic',T=5,ratio=ratio)
    print(sc,(['%.2e'%e for e in E],['%.2f'%p for p in orders(E)[1:]]) if E else 'X','drift',['%.1e'%r.get('drift',float('nan')) for r in res] )
    out[f'ms_{sc}']=(E,orders(E)[1:] if E else None,[r['drift'] for r in res])
json.dump(out,open(os.path.join(ROOT,'results','exp3.json'),'w'),default=lambda o: float(o) if np.isscalar(o) else str(o))
