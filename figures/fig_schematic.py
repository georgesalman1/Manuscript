import os, numpy as np, matplotlib; matplotlib.use('Agg')
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.spatial import Voronoi
plt.rcParams.update({'font.size':9,'font.family':'serif','mathtext.fontset':'cm'})
fig,(ax1,ax2)=plt.subplots(2,1,figsize=(6.6,4.4),gridspec_kw={'height_ratios':[1.5,1]})
n,ny,dx,dy=8,3,1.0,0.75
I,J=np.meshgrid(np.arange(-1,n+3),np.arange(-2,ny+2),indexing='ij')
pts=np.c_[I.ravel()*dx,(J.ravel()+0.5)*dy]
vor=Voronoi(pts)
inside=lambda p:(0.5<=p[0]<=n+0.5) and (0<=p[1]<=ny*dy)
for (p,q),rv in zip(vor.ridge_points,vor.ridge_vertices):
    if -1 in rv: continue
    a,b=pts[p],pts[q]
    if inside(a) and inside(b):
        V=vor.vertices[rv]
        if np.linalg.norm(V[1]-V[0])<1e-9: continue
        ax1.plot(V[:,0],V[:,1],'k-',lw=0.8)
ax1.add_patch(Rectangle((0.5,0),n,ny*dy,fill=False,lw=1.6))
cols={3:'#4C72B0',4:'#C44E52',5:'#55A868'}
for i in range(1,n+1):
    for j in range(ny):
        c=cols.get(i,'0.4')
        ax1.plot(i*dx,(j+.5)*dy,'o',ms=4,color=c)
for i,lab in [(3,'$l$'),(4,'$m$'),(5,'$r$')]:
    ax1.add_patch(Rectangle((i-0.5,dy),1,dy,color=cols[i],alpha=0.18,lw=0))
    ax1.text(i,1.5*dy+0.13,lab,ha='center',va='bottom',fontsize=10,color=cols[i])
# normals on cell m
y0=1.5*dy
ax1.annotate('',xy=(4.5+0.38,y0),xytext=(4.5,y0),arrowprops=dict(arrowstyle='->',lw=1.2))
ax1.annotate('',xy=(4.5-1-0.38,y0),xytext=(3.5,y0),arrowprops=dict(arrowstyle='->',lw=1.2))
ax1.annotate('',xy=(4,2*dy+0.3),xytext=(4,2*dy),arrowprops=dict(arrowstyle='->',lw=1.0,color='0.5'))
ax1.text(4.9,y0-0.17,r'$\mathbf{n}_{mr}$',fontsize=8); ax1.text(3.0,y0-0.17,r'$\mathbf{n}_{ml}$',fontsize=8)
ax1.text(4.08,2*dy+0.16,r'$\mathbf{u}\cdot\mathbf{n}_{ij}=0$ (transverse faces)',fontsize=7,color='0.4')
ax1.text(0.6,ny*dy+0.12,'solid wall: zero normal flux, hydrostatic wall pressure',fontsize=7.5)
ax1.text(0.55,-0.32,'Dirichlet column $i=1$',fontsize=7.5); ax1.text(n-1.6,-0.32,'Dirichlet column $i=n$',fontsize=7.5)
ax1.set_xlim(0.2,n+0.8); ax1.set_ylim(-0.5,ny*dy+0.45); ax1.set_aspect('equal'); ax1.axis('off')
ax1.set_title('(a) Quadrilateral Voronoi region: generators, faces $A_{ij}$, unit normals $\\mathbf{n}_{ij}$',fontsize=9,loc='left')
# compressed
for i in range(1,n+1):
    c=cols.get(i,'0.4')
    ax2.add_patch(Rectangle((i-0.5,0),1,0.6,fill=i in cols,color=c if i in cols else None,alpha=0.18 if i in cols else 1,lw=0))
    ax2.add_patch(Rectangle((i-0.5,0),1,0.6,fill=False,lw=0.8))
    ax2.plot(i,0.3,'o',ms=4,color=c)
    ax2.text(i,-0.18,str(i) if i not in (n,) else '$n$',ha='center',fontsize=8)
for i,lab in [(3,'$l$'),(4,'$m$'),(5,'$r$')]:
    ax2.text(i,0.7,lab,ha='center',fontsize=10,color=cols[i])
def brace(x0,x1,y,txt):
    ax2.annotate('',xy=(x0,y),xytext=(x1,y),arrowprops=dict(arrowstyle='|-|',lw=0.8,mutation_scale=3))
    ax2.text((x0+x1)/2,y-0.2,txt,ha='center',fontsize=8)
brace(1.5,n-0.5,-0.42,'$m=2{:}n-1$  (updated)')
ax2.text(0.55,-0.85,r'$a=1{:}n,\;\; l=1{:}n-2,\;\; m=2{:}n-1,\;\; r=3{:}n$;  $\;\Delta x = V_i/A_{i,i+1}$',fontsize=8.5)
ax2.set_xlim(0.2,n+0.8); ax2.set_ylim(-1.0,0.95); ax2.set_aspect('equal'); ax2.axis('off')
ax2.set_title('(b) Compressed one-dimensional topology and index sets',fontsize=9,loc='left')
fig.tight_layout(); fig.savefig(os.path.join(ROOT,'results','figs','compression_schematic.pdf')); fig.savefig(os.path.join(ROOT,'results','figs','compression_schematic.png'),dpi=110)
