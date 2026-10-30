import os, sys, json, numpy as np, matplotlib; matplotlib.use('Agg'); ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,os.path.join(ROOT,'src'))
import os, matplotlib.pyplot as plt
from mmg import *
plt.rcParams.update({'font.size':9,'font.family':'serif','mathtext.fontset':'cm','axes.grid':True,'grid.alpha':0.3})
fig,ax=plt.subplots(1,3,figsize=(7.4,2.7))
# (a) k sweep PA (from sweep log)
ks={1.05:[4.38e-03,2.36e-03,1.24e-03,6.37e-04],1.25:[2.44e-03,1.27e-03,6.49e-04,3.28e-04],1.5:[1.34e-03,6.84e-04,3.45e-04,1.74e-04],
    2:[5.55e-04,2.77e-04,1.38e-04,6.89e-05],3:[1.96e-04,9.42e-05,4.61e-05,2.28e-05],5:[7.15e-05,3.26e-05,1.55e-05,7.54e-06]}
dx=np.array([1/50,1/100,1/200,1/400]); cm=plt.cm.viridis
for i,(k,E) in enumerate(ks.items()):
    ax[0].loglog(dx,E,'o-',ms=3.5,color=cm(i/5),label=f'$k={k}$')
ax[0].loglog(dx,0.1*dx,'k--',lw=0.8)
ax[0].set_xlabel(r'$\Delta x$'); ax[0].set_ylabel(r'steady $L_2$ error'); ax[0].legend(fontsize=6.5,ncol=2)
ax[0].set_title('(a) mean depth (PA, MMG)',fontsize=9,loc='left')
ax[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[0].set_xticks([2.5e-3,1e-2]); ax[0].set_xticklabels(['0.0025','0.01'])
# (b) multiscale
M=multiscale_profile()
o=run(200,'pa',M,bc='periodic',T=5,ratio=0.05)
xf=np.linspace(0,1,1000)
ax[1].plot(xf,-1/M.s(xf),'k-',lw=0.9,label='manufactured')
ax[1].plot(o['x'][::2],o['u'][::2],'o',ms=2.5,mfc='none',color='#4C72B0',label='PA, $n=200$')
ax[1].set_xlabel('$x$'); ax[1].set_ylabel('$u$'); ax[1].legend(fontsize=6.5,loc='lower right')
ax[1].set_title('(b) multiscale profile',fontsize=9,loc='left')
# (c) long time
d=json.load(open(os.path.join(ROOT,'results','exp2.json')))
C={'pa':'#4C72B0','cu':'#C44E52','mms':'#55A868'}; L={'pa':'PA (MMG)','cu':'CU (MMG)','mms':'MMS'}
for sc in C:
    h=np.array(d[f'long_{sc}']); ax[2].loglog(h[:,0],h[:,1],'-',color=C[sc],label=L[sc])
ax[2].set_xlabel('$t$'); ax[2].set_ylabel(r'$L_2$ error, $n=100$'); ax[2].legend(fontsize=6.5)
ax[2].set_title('(c) long integration',fontsize=9,loc='left')
fig.tight_layout(); fig.savefig(os.path.join(ROOT,'results','figs','robustness.pdf')); fig.savefig(os.path.join(ROOT,'results','figs','robustness.png'),dpi=90)
