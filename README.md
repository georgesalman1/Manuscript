# MMG/MMS verification scripts

Supplementary code for: G. Salman, C. Liosis, S. Peppa, I. Sarris, "Manufactured Geometry vs. Manufactured Forcing: Verifying a Fluid Dynamic Voronoi Model", *Symmetry* (MDPI), manuscript symmetry-4631884.

This is a minimal, standalone Python implementation of the one-dimensional verification tests (MMG and MMS) described in the paper. It is not the Wavedyne code itself; it implements the same discrete stencils on the compressed one-dimensional topology, plus a generic two-dimensional Voronoi face-flux check (Proposition 1).

It reproduces every table and figure of the paper.

Requirements: Python 3 (tested with Python 3.13, NumPy 2.5, SciPy 1.18, Matplotlib 3.11). Install with `pip install -r requirements.txt`.

| Script | Produces |
|---|---|
| `mmg.py` | Solver library: MMG (product-average / convective upwind), MMS (continuous / discrete forcing), mirrored and sign-adaptive stencils, periodic BCs, fault injection, discrete residual |
| `exp1.py` | Baseline Table 1, discrete MMS forcing, discrete residual, steady state vs T=2, volume budget |
| `exp2.py` | Perturbed initial states, long integration T=200 (writes exp2.json, used by fig_robust.py) |
| `exp3.py` | k sweep (Table 3), rightward flow, periodic k=3 and multiscale profiles (Table 4) |
| `exp4.py` | Fault injection (Table 5), sign-changing MMS fields |
| `exp5.py` | Stability limit vs k |
| `exp6.py` | Multiscale MMS time-step scan |
| `exp7.py` | k=5 at dt/dx=0.02, multiscale PA, periodic momentum conservation |
| `voronoi2d.py` | 2D Voronoi strip vs compressed 1D equivalence (Proposition 1) |
| `symmetry_check.py` | Section 2.3.5: skew-symmetric/symmetric decomposition, odd-even damping, momentum conservation |
| `exp8.py` | Periodic product-average case run to true steady state (Table 4) |
| `fig_*.py` | Figures 3, 4, 5, 6, 7, 8 |

Run from this folder, e.g. `python3 exp1.py`. The k sweep (exp3.py) takes about 10–15 minutes.

Note on periodic tests: the T=40 periodic results printed by exp3.py are transient. The steady-state values in Table 4 come from `exp8.py N` (N = 50, 100, 200), which integrates the periodic product-average case to a residual below 1e-11 (t ≈ 3.2e3, 6.2e3 and 1.17e4; n = 200 takes about 15 minutes).
