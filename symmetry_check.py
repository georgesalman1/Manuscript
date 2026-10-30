"""Section 2.3.5: numerical check of the skew-symmetric / symmetric decomposition (Eq. 18),
odd-even damping rate 4 nu/dx^2, and momentum conservation of the product-average stencil."""
import numpy as np
n=16; ub=-0.7; dx=0.1; eps=1e-7
I=np.arange(n); R=np.roll(I,-1); L=np.roll(I,1)
CU=lambda u: u*(u[R]-u)/dx
PA=lambda u: 0.25*((u[R]+u)*u[R]-(u+u[L])*u)/dx
def jac(f):
    J=np.zeros((n,n)); u0=np.full(n,ub)
    for j in range(n):
        e=np.zeros(n); e[j]=eps; J[:,j]=(f(u0+e)-f(u0-e))/(2*eps)
    return J
Dc=np.zeros((n,n)); D2=np.zeros((n,n))
for i in range(n):
    Dc[i,R[i]]+=1/(2*dx); Dc[i,L[i]]-=1/(2*dx); D2[i,R[i]]+=1/dx**2; D2[i,L[i]]+=1/dx**2; D2[i,i]-=2/dx**2
for name,f,nu in [('CU',CU,abs(ub)*dx/2),('PA',PA,abs(ub)*dx/4)]:
    J=jac(f); v=(-1.0)**I
    print(name,'max|J-(ub*Dc-nu*D2)| = %.1e'%abs(J-(ub*Dc-nu*D2)).max(),
          ' odd-even damping %.4f (predicted %.4f)'%(v@J@v/(v@v),4*nu/dx**2))
u=np.random.default_rng(0).random(n)+0.5
print('sum over periodic grid: PA %.1e, CU %.1e'%(PA(u).sum(),CU(u).sum()))
