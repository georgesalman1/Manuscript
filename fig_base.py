import sys, json, numpy as np, matplotlib; matplotlib.use('Agg'); sys.path.insert(0,'.')
import matplotlib.pyplot as plt
from mmg import *
plt.rcParams.update({'font.size':9,'font.family':'serif','mathtext.fontset':'cm','axes.grid':True,'grid.alpha':0.3})
C={'pa':'#4C72B0','cu':'#C44E52','mms':'#55A868'}; LAB={'pa':'Product-average (MMG)','cu':'Convective upwind (MMG)','mms':'MMS scalar advection'}
MK={'pa':'o','cu':'s','mms':'^'}
P=sine_profile(2.0); N=[50,100,200,400]
R={sc:[run(n,sc,P) for n in N] for sc in C}
# --- convergence
fig,ax=plt.subplots(figsize=(4.2,3.3))
dx=np.array([1/n for n in N])
for sc in C:
    ax.loglog(dx,[r['L2'] for r in R[sc]],MK[sc]+'-',color=C[sc],label=LAB[sc],mfc='none' if sc=='pa' else C[sc],ms=5)
ax.loglog(dx,0.03*dx,'k--',lw=0.8,label=r'$O(\Delta x)$')
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xticks([2.5e-3,5e-3,1e-2,2e-2]); ax.set_xticklabels(['0.0025','0.005','0.01','0.02']); ax.set_xlabel(r'$\Delta x$'); ax.set_ylabel(r'$L_2$ error in $u$ at $T_{final}=2$'); ax.legend(fontsize=7.5)
fig.tight_layout(); fig.savefig('figs/convergence.pdf')
# --- velocity + error
fig,(a1,a2)=plt.subplots(2,1,figsize=(5.2,4.6),sharex=True,gridspec_kw={'height_ratios':[1.4,1]})
r=R['pa'][1]; xf=np.linspace(0.005,1,400)
a1.plot(xf,-1/(2+np.sin(xf)),'k-',lw=1.2,label=r'manufactured $u=-s^{-1}$')
for sc in C:
    o=R[sc][1]; a1.plot(o['x'][::6],o['u'][::6],MK[sc],color=C[sc],mfc='none',ms=5,label=LAB[sc])
    a2.semilogy(o['x'][1:-1],abs(o['u']-o['ue'])[1:-1],'-',color=C[sc],label=LAB[sc])
a1.set_ylabel('$u$'); a1.legend(fontsize=7.5)
a2.set_ylabel(r'$|u_{num}-u_{exact}|$'); a2.set_xlabel('$x$'); a2.set_ylim(1e-5,2e-3)
fig.tight_layout(); fig.savefig('figs/velocity_error.pdf')
# --- domain
fig,ax=plt.subplots(figsize=(5.6,3.2))
o=R['pa'][1]; x=o['x']; b=o['b']
ax.fill_between(x,-b,-0.5*(1/(2+np.sin(x)))**2,color='#4C72B0',alpha=0.15,lw=0,label='MMG fluid layer')
ax.plot(x,-b,'k-',lw=1.3,label=r'MMG bed $z=-b(x)$')
ax.plot(xf,-0.5/(2+np.sin(xf))**2,'-',color='#4C72B0',lw=1,label=r'exact MMG surface $z=e=-u^2/2$')
for sc in ['pa','cu']:
    oo=R[sc][1]; ax.plot(oo['x'][::8],(oo['h']-oo['b'])[::8],MK[sc],color=C[sc],mfc='none',ms=5,label='computed surface, '+LAB[sc].split(' (')[0].lower())
ax.plot([0,1],[-2,-2],'--',color=C['mms'],lw=1.2,label=r'MMS bed $z=-k$')
ax.plot([0,1],[0,0],':',color=C['mms'],lw=1.4,label=r'MMS surface $z=0$')
ax.set_xlabel('$x$'); ax.set_ylabel('elevation $z$'); ax.set_ylim(-3.1,0.35); ax.legend(fontsize=6.8,loc='center',ncol=2,framealpha=0.95)
fig.tight_layout(); fig.savefig('figs/domain.pdf')
# --- residual decomposition + dt sweep
fig,(a1,a2)=plt.subplots(1,2,figsize=(7.0,2.9))
NN=[50,100,200,400,800]
D={sc:[discrete_residual(n,sc,P) for n in NN] for sc in ['pa','cu']}
d=np.array([1/n for n in NN])
a1.loglog(d,[q['adv_trunc'] for q in D['pa']],'o-',color=C['pa'],mfc='none',label='advective truncation (PA)')
a1.loglog(d,[q['adv_trunc'] for q in D['cu']],'s-',color=C['cu'],label='advective truncation (CU)')
a1.loglog(d,[q['grav_trunc'] for q in D['pa']],'d-',color='0.4',label='hydrostatic truncation')
a1.loglog(d,[q['Rh'] for q in D['pa']],'v-',color='#8172B2',label=r'continuity residual $R_h$')
a1.loglog(d,0.08*d,'k--',lw=0.8,label=r'$O(\Delta x)$'); a1.loglog(d,0.04*d**2,'k:',lw=0.8,label=r'$O(\Delta x^2)$')
a1.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); a1.set_xticks([1.25e-3,5e-3,2e-2]); a1.set_xticklabels(['0.00125','0.005','0.02']); a1.set_xlabel(r'$\Delta x$'); a1.set_ylabel(r'$L_\infty$ residual at exact state'); a1.legend(fontsize=6.5); a1.set_title('(a) discrete steady residual',fontsize=9,loc='left')
ratios=[0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18]
lim={'pa':0.1361,'cu':0.1802}
for sc in ['pa','cu']:
    E=[]
    for q in ratios:
        o=run(100,sc,P,ratio=q,T=20); E.append(np.nan if o['blowup'] else o['L2'])
    a2.semilogy(ratios,E,MK[sc]+'-',color=C[sc],mfc='none' if sc=='pa' else C[sc],label=LAB[sc])
    a2.axvline(lim[sc],color=C[sc],ls='--',lw=0.8)
a2.set_xlabel(r'$\Delta t/\Delta x$'); a2.set_ylabel(r'$L_2$ error, $n=100$, $T=20$'); a2.set_ylim(1e-4,1e-3)
a2.legend(fontsize=6.8); a2.set_title(r'(b) error vs. time step',fontsize=9,loc='left')
fig.tight_layout(); fig.savefig('figs/residual_dt.pdf')
print('ok')
