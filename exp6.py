import sys, json, numpy as np; sys.path.insert(0,'.')
from mmg import *
out={}
M=multiscale_profile(); N2=[50,100,200,400,800,1600]
print("== MMS multiscale, ratio scan n=200")
for ratio in [0.5,0.2,0.1,0.05]:
    o=run(200,'mms',M,bc='periodic',T=5,ratio=ratio); print(ratio,'blowup' if o['blowup'] else '%.3e'%o['L2'], o.get('t'))
