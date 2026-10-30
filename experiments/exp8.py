import os, sys, numpy as np; ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,os.path.join(ROOT,'src'))
from mmg import *
P=sine_profile(3.0); n=int(sys.argv[1])
o=run(n,'pa',P,bc='periodic',L=2*np.pi,to_steady=True,steady_tol=1e-11,ratio=0.04,T=1e9)
print(n,'L2 %.4e t_ss %.0f flux %.6f drift %.1e'%(o['L2'],o['t'],np.mean(o['u']*o['h']),o['drift']),flush=True)
