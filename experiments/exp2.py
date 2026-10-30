import os, sys, json, numpy as np; ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,os.path.join(ROOT,'src'))
from mmg import *
P=sine_profile(2.0); N=[50,100,200,400]
out={}
print("== perturbed initial state (10% of |u|), run to discrete steady state")
for sc,ratio in [('pa',0.1),('cu',0.1),('mms',0.5)]:
    for forc in (['continuous','discrete'] if sc=='mms' else ['continuous']):
        row=[]
        for n in N:
            o=run(n,sc,P,perturb=0.1,to_steady=True,ratio=ratio,forcing=forc,steady_tol=1e-12)
            o0=run(n,sc,P,to_steady=True,ratio=ratio,forcing=forc,steady_tol=1e-12)
            row.append((n,o['L2'],o0['L2'],abs(o['u']-o0['u']).max(),o['t']))
        print(sc,forc,[('%d'%a,'%.3e'%b,'%.3e'%c,'%.1e'%d,'%.1f'%t) for a,b,c,d,t in row])
        out[f'pert_{sc}_{forc}']=row
print("== long time T=200, n=100, ratio 0.02")
rec=sorted(set(np.unique(np.logspace(0,6,61).astype(int))))
for sc in ['pa','cu','mms']:
    o=run(100,sc,P,T=200,record=rec)
    print(sc,'L2(T=200) %.4e drift %.2e'%(o['L2'],o['drift']), 'hist first/last',o['hist'][0],o['hist'][-1])
    out[f'long_{sc}']=o['hist']; out[f'long_{sc}_L2']=o['L2']; out[f'long_{sc}_drift']=o['drift']
json.dump(out,open(os.path.join(ROOT,'results','exp2.json'),'w'),default=float)
