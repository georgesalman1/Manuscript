import numpy as np, json, sys
sys.path.insert(0,'.')
from mmg import *
P=sine_profile(2.0); N=[50,100,200,400]
print("== baseline")
for sc in ['pa','cu','mms']:
    E=[run(n,sc,P)['L2'] for n in N]; print(sc,['%.3e'%e for e in E],['%.3f'%p for p in orders(E)[1:]])
print("== discrete MMS forcing")
for n in N:
    o=run(n,'mms',P,forcing='discrete'); print(n,'%.2e %.2e'%(o['L2'],o['Linf']))
print("== discrete residual at exact state")
for sc in ['pa','cu']:
    for n in N:
        d=discrete_residual(n,sc,P); print(sc,n,{k:'%.2e'%v for k,v in d.items()})
print("== steady state vs T=2")
for sc in ['pa','cu','mms']:
    for n in [50,100]:
        a=run(n,sc,P); s=run(n,sc,P,to_steady=True,ratio=0.2 if sc!='mms' else 0.2)
        print(sc,n,'T2 %.4e  steady %.4e  t_ss %.1f'%(a['L2'],s.get('L2',np.nan),s.get('t',np.nan)), s['blowup'])
print("== volume budget")
for sc in ['pa','cu']:
    o=run(100,sc,P); print(sc,'vol change %.6e  pred from boundary fluxes %.6e'%(o['vol_change'],o['bflux_pred']))
