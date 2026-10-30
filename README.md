# Manufactured Geometry vs. Manufactured Forcing: MMG/MMS verification code

Supplementary code for:

> G. Salman, C. Liosis, S. Peppa, I. Sarris, "Manufactured Geometry vs. Manufactured Forcing: Verifying a Fluid Dynamic Voronoi Model", *Symmetry* (MDPI), 2026. Under review, manuscript symmetry-4631884.
>
> **Publication:** link will be added once the paper is published (DOI pending).

This is a minimal, standalone Python implementation of the one-dimensional verification tests (MMG and MMS) described in the paper. It is not the Wavedyne code itself; it implements the same discrete stencils on the compressed one-dimensional topology, plus a generic two-dimensional Voronoi face-flux check (Proposition 1).

It reproduces every table and figure of the paper.

## Repository layout

```
src/            solver library (mmg.py)
experiments/    scripts that produce the tables and numerical checks
figures/        scripts that produce the paper figures
results/        generated output: exp*.json and figs/ (not tracked by git)
```

## Installation

Python 3 (tested with Python 3.13, NumPy 2.5, SciPy 1.18, Matplotlib 3.11).

```
pip install -r requirements.txt
```

## Usage

Scripts can be run from any directory, e.g.

```
python3 experiments/exp1.py
python3 figures/fig_base.py
```

JSON output goes to `results/`, figures to `results/figs/`. Run `experiments/exp2.py` before `figures/fig_robust.py`, which reads `results/exp2.json`.

| Script | Produces |
|---|---|
| `src/mmg.py` | Solver library: MMG (product-average / convective upwind), MMS (continuous / discrete forcing), mirrored and sign-adaptive stencils, periodic BCs, fault injection, discrete residual |
| `experiments/exp1.py` | Baseline Table 1, discrete MMS forcing, discrete residual, steady state vs T=2, volume budget |
| `experiments/exp2.py` | Perturbed initial states, long integration T=200 (writes `results/exp2.json`, used by `fig_robust.py`) |
| `experiments/exp3.py` | k sweep (Table 3), rightward flow, periodic k=3 and multiscale profiles (Table 4) |
| `experiments/exp4.py` | Fault injection (Table 5), sign-changing MMS fields |
| `experiments/exp5.py` | Stability limit vs k |
| `experiments/exp6.py` | Multiscale MMS time-step scan |
| `experiments/exp7.py` | k=5 at dt/dx=0.02, multiscale PA, periodic momentum conservation |
| `experiments/exp8.py` | Periodic product-average case run to true steady state (Table 4) |
| `experiments/voronoi2d.py` | 2D Voronoi strip vs compressed 1D equivalence (Proposition 1) |
| `experiments/symmetry_check.py` | Section 2.3.5: skew-symmetric/symmetric decomposition, odd-even damping, momentum conservation |
| `figures/fig_*.py` | Figures 3, 4, 5, 6, 7, 8 |

The k sweep (`exp3.py`) takes about 10–15 minutes.

Note on periodic tests: the T=40 periodic results printed by `exp3.py` are transient. The steady-state values in Table 4 come from `experiments/exp8.py N` (N = 50, 100, 200), which integrates the periodic product-average case to a residual below 1e-11 (t ≈ 3.2e3, 6.2e3 and 1.17e4; n = 200 takes about 15 minutes).

## How to cite

If you use this code, please cite the paper above. Citation metadata is in [`CITATION.cff`](CITATION.cff) (GitHub shows it under "Cite this repository").

## License

MIT, see [`LICENSE`](LICENSE).
