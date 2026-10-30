"""Generic 2D Voronoi face-flux discretization of the MMG system, evaluated on a strip of
quadrilateral Voronoi cells, compared with the compressed 1D update.
All geometry (areas, face lengths, normals, neighbor pairs) is taken from scipy's Voronoi
construction, not assumed."""
import os, sys
import numpy as np
from scipy.spatial import Voronoi
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,os.path.join(ROOT,'src'))
from mmg import sine_profile, run

g = 1.0


def build(n, ny, dx, dy, jitter=0.0, seed=0):
    rng = np.random.default_rng(seed)
    # interior generators + two ghost layers on every side to close the regions
    I, J = np.meshgrid(np.arange(-1, n + 3), np.arange(-2, ny + 2), indexing='ij')
    pts = np.c_[I.ravel() * dx, (J.ravel() + 0.5) * dy]
    interior = ((I >= 1) & (I <= n) & (J >= 0) & (J < ny)).ravel()
    pts = pts + jitter * dx * rng.uniform(-1, 1, pts.shape) * interior[:, None]
    vor = Voronoi(pts)
    idmap = -np.ones(len(pts), int)
    idmap[np.where(interior)[0]] = np.arange(interior.sum())
    cells = pts[interior]
    col = I.ravel()[interior]
    # areas by shoelace on the Voronoi regions
    area = np.empty(len(cells))
    for k, p in enumerate(np.where(interior)[0]):
        reg = vor.regions[vor.point_region[p]]
        V = vor.vertices[reg]
        c = V.mean(0)
        ang = np.arctan2(V[:, 1] - c[1], V[:, 0] - c[0])
        V = V[np.argsort(ang)]
        area[k] = 0.5 * abs(np.dot(V[:, 0], np.roll(V[:, 1], -1)) - np.dot(V[:, 1], np.roll(V[:, 0], -1)))
    faces = []
    for (p, q), rv in zip(vor.ridge_points, vor.ridge_vertices):
        a, b = idmap[p], idmap[q]
        if a < 0 or b < 0 or -1 in rv:
            continue                      # wall / ghost faces carry no flux
        V = vor.vertices[rv]
        Lf = np.linalg.norm(V[1] - V[0])
        if Lf < 1e-12:
            continue                      # degenerate (corner-touching) ridges
        nrm = (cells[b] - cells[a]) / np.linalg.norm(cells[b] - cells[a])
        faces.append((a, b, Lf, nrm[0], nrm[1]))
    F = np.array(faces)
    return cells, col, area, F


def run2d(n, ny, scheme, T=2.0, ratio=0.02, k=2.0):
    dx = 1.0 / n; dy = 0.7 * dx; dt = ratio * dx; Nt = int(round(T / dt))
    cells, col, A, F = build(n, ny, dx, dy)
    a, b, Lf, nx, ny_ = F[:, 0].astype(int), F[:, 1].astype(int), F[:, 2], F[:, 3], F[:, 4]
    x = cells[:, 0]
    s = k + np.sin(x); ue = -1 / s; bed = s + 0.5 * ue**2
    u = ue.copy(); v = np.zeros_like(u); h = s.copy()
    active = (col > 1) & (col < n)        # Dirichlet-anchored first and last columns
    # wall faces: closure defect D_a = -sum_f L n (missing area vector of solid-wall faces);
    # hydrostatic wall pressure e_wall = e_a contributes e_a * D_a to the Green-Gauss sum
    Dx = -(np.bincount(a, Lf * nx, len(u)) - np.bincount(b, Lf * nx, len(u)))
    Dy = -(np.bincount(a, Lf * ny_, len(u)) - np.bincount(b, Lf * ny_, len(u)))
    for _ in range(Nt):
        e = h - bed
        Un = 0.5 * ((u[a] + u[b]) * nx + (v[a] + v[b]) * ny_)     # face-normal velocity a->b
        up = Un < 0                                                # inflow into a -> upwind is b
        # continuity: outward flux from a
        fh = Lf * Un * np.where(up, h[b], h[a])
        divh = np.bincount(a, fh, len(u)) - np.bincount(b, fh, len(u))
        # hydrostatic Green-Gauss gradient
        ef = 0.5 * (e[a] + e[b]) * Lf
        gx = np.bincount(a, ef * nx, len(u)) - np.bincount(b, ef * nx, len(u))
        gy = np.bincount(a, ef * ny_, len(u)) - np.bincount(b, ef * ny_, len(u)) + e * Dy
        gx = gx + e * Dx * active
        if scheme == 'pa':
            fu = Lf * Un * 0.5 * np.where(up, u[b], u[a])
            fv = Lf * Un * 0.5 * np.where(up, v[b], v[a])
            advu = np.bincount(a, fu, len(u)) - np.bincount(b, fu, len(u))
            advv = np.bincount(a, fv, len(u)) - np.bincount(b, fv, len(u))
        else:   # convective upwind: sum over inflow faces of (u_i.n) (phi_j - phi_i) L
            una = u[a] * nx + v[a] * ny_            # cell-a velocity along a->b
            unb = -(u[b] * nx + v[b] * ny_)         # cell-b velocity along b->a
            wa = np.where(una < 0, Lf * una, 0.0)   # face is inflow for a
            wb = np.where(unb < 0, Lf * unb, 0.0)
            advu = np.bincount(a, wa * (u[b] - u[a]), len(u)) + np.bincount(b, wb * (u[a] - u[b]), len(u))
            advv = np.bincount(a, wa * (v[b] - v[a]), len(u)) + np.bincount(b, wb * (v[a] - v[b]), len(u))
        hn = h - dt * divh / A
        un = u - dt * (advu + g * gx) / A
        vn = v - dt * (advv + g * gy) / A
        h = np.where(active, hn, h); u = np.where(active, un, u); v = np.where(active, vn, v)
    return dict(x=x, col=col, u=u, v=v, h=h, A=A, nfaces=len(F), dx=dx, dy=dy)


if __name__ == '__main__':
    P = sine_profile(2.0)
    for n in [50, 100]:
        for sc in ['pa', 'cu']:
            r2 = run2d(n, 4, sc)
            r1 = run(n, sc, P)
            du = abs(r2['u'] - r1['u'][r2['col'] - 1]).max()
            dh = abs(r2['h'] - r1['h'][r2['col'] - 1]).max()
            print(n, sc, 'faces', r2['nfaces'], 'max|u2D-u1D| %.2e  max|h2D-h1D| %.2e  max|v| %.2e  area/dxdy-1 %.1e'
                  % (du, dh, abs(r2['v']).max(), abs(r2['A'] / (r2['dx'] * r2['dy']) - 1).max()))
