"""General 1D solver reproducing the Wavedyne verification runs (MMG and MMS)
and the additional configurations requested in the review."""
import numpy as np

g = 1.0


class Profile:
    """Manufactured depth profile s(x) with derivative ds(x)."""
    def __init__(self, s, ds, name=""):
        self.s, self.ds, self.name = s, ds, name


def sine_profile(k):
    return Profile(lambda x: k + np.sin(x), lambda x: np.cos(x), f"k+sin(x), k={k}")


def multiscale_profile(k=2.0, a=(0.5, 0.25, 0.1), w=(2, 8, 32)):
    def s(x):
        return k + sum(ai * np.sin(wi * np.pi * x) for ai, wi in zip(a, w))
    def ds(x):
        return sum(ai * wi * np.pi * np.cos(wi * np.pi * x) for ai, wi in zip(a, w))
    return Profile(s, ds, "multiscale")


def run(n, scheme, prof, C=-1.0, direction="fwd", bc="dirichlet", forcing="continuous",
        T=2.0, ratio=0.02, L=1.0, mutation=None, record=None, uexact=None, duexact=None,
        to_steady=False, steady_tol=1e-13, max_steps=None, perturb=0.0):
    """scheme: 'pa' product-average MMG, 'cu' convective upwind MMG, 'mms' scalar advection.
    direction: 'fwd' (original forward-biased, upwind for u<0), 'bwd' (mirrored, upwind for u>0),
               'adaptive' (sign-switching upwind, MMS only).
    forcing (MMS): 'continuous' or 'discrete'.
    mutation: None or one of 'grav2' (missing 1/2 in hydrostatic term), 'bedsign' (e=h+b),
              'adv2dx' (advection divided by 2dx), 'hcentral' (depth flux uses h at own cell)."""
    dx = L / n
    dt = ratio * dx
    Nt = int(round(T / dt)) if max_steps is None else max_steps
    x = np.arange(1, n + 1) * dx
    s = prof.s(x)
    if uexact is None:
        ue = C / s
        due = -C * prof.ds(x) / s**2
    else:
        ue, due = uexact(x), duexact(x)
    b = s + 0.5 * ue**2 / g          # manufactured bed, b = h + u^2/(2g)
    u = ue.copy()
    if perturb:
        u[1:-1] += perturb * np.sin(np.pi * x[1:-1] / L) * np.abs(ue[1:-1])
    if bc == "periodic":
        idx = np.arange(n)
        m, l, r = idx, np.roll(idx, 1), np.roll(idx, -1)
    else:
        m, l, r = np.arange(1, n - 1), np.arange(0, n - 2), np.arange(2, n)
    if scheme == "mms":
        h = np.full(n, 2.0)
        if forcing == "continuous":
            F = ue * due
        else:
            if direction == "fwd":
                F = np.zeros(n); F[m] = ue[m] * (ue[r] - ue[m]) / dx
            elif direction == "bwd":
                F = np.zeros(n); F[m] = ue[m] * (ue[m] - ue[l]) / dx
            else:
                F = np.zeros(n)
                F[m] = np.where(ue[m] < 0, ue[m] * (ue[r] - ue[m]), ue[m] * (ue[m] - ue[l])) / dx
    else:
        h = s.copy()
    sign_b = -1.0 if mutation == "bedsign" else 1.0
    gfac = 1.0 if mutation == "grav2" else 0.5
    advfac = 0.5 if mutation == "adv2dx" else 1.0
    V0 = h.sum()
    drift = 0.0
    bflux_int = 0.0
    hist = []
    rec_set = set(record) if record is not None else set()

    def adv_of(u):
        if direction == "fwd":
            if scheme == "pa":
                return 0.25 * ((u[r] + u[m]) * u[r] - (u[m] + u[l]) * u[m]) / dx
            return u[m] * (u[r] - u[m]) / dx
        if direction == "bwd":
            if scheme == "pa":
                return 0.25 * ((u[r] + u[m]) * u[m] - (u[m] + u[l]) * u[l]) / dx
            return u[m] * (u[m] - u[l]) / dx
        # adaptive
        return np.where(u[m] < 0, u[m] * (u[r] - u[m]), u[m] * (u[m] - u[l])) / dx

    def hflux(u, h):
        # face fluxes F_{m+1/2}, F_{m-1/2}
        if mutation == "hcentral":
            return 0.5 * (u[r] + u[m]) * 0.5 * (h[r] + h[m]), 0.5 * (u[m] + u[l]) * 0.5 * (h[m] + h[l])
        if direction == "fwd":
            return 0.5 * (u[r] + u[m]) * h[r], 0.5 * (u[m] + u[l]) * h[m]
        return 0.5 * (u[r] + u[m]) * h[m], 0.5 * (u[m] + u[l]) * h[l]

    step = 0
    while True:
        if scheme == "mms":
            res_u = advfac * adv_of(u) - F[m]
            un = u.copy(); un[m] = u[m] - dt * res_u
            res_h = 0.0
            u = un
        else:
            e = h - sign_b * b
            Fp, Fm = hflux(u, h)
            res_h = (Fp - Fm) / dx
            res_u = advfac * adv_of(u) + g * gfac * (e[r] - e[l]) / dx
            if bc != "periodic":
                bflux_int += dt * (Fp[-1] - Fm[0])
            hn = h.copy(); hn[m] = h[m] - dt * res_h
            un = u.copy(); un[m] = u[m] - dt * res_u
            u, h = un, hn
            drift = max(drift, abs(h.sum() - V0) / V0)
        step += 1
        if not np.all(np.isfinite(u)) or np.abs(u).max() > 1e6:
            return dict(blowup=True, step=step, t=step * dt)
        if step in rec_set:
            ee = u - ue
            hist.append((step * dt, np.sqrt(np.mean(ee**2))))
        if to_steady:
            rr = max(np.abs(res_u).max(), np.abs(res_h).max() if scheme != "mms" else 0)
            if rr < steady_tol or step >= 50_000_000:
                break
        elif step >= Nt:
            break
    eps = u - ue
    out = dict(blowup=False, n=n, dx=dx, L1=np.mean(abs(eps)), L2=np.sqrt(np.mean(eps**2)),
               Linf=abs(eps).max(), drift=drift, x=x, u=u, ue=ue, h=h, b=b, steps=step,
               t=step * dt, hist=hist)
    if scheme != "mms":
        out["vol_change"] = (h.sum() - V0)
        out["bflux_int"] = -bflux_int / dx * dx  # integrated net boundary flux (per dx units)
        # sum h changes by -dt/dx * sum(Fp-Fm) = -dt/dx*(F_{n-1/2}-F_{3/2})
        out["bflux_pred"] = -bflux_int / dx
    cs = np.abs(ue) + np.sqrt(g * s)
    out["cfl"] = cs.max() * ratio
    out["Fr_max"] = (np.abs(ue) / np.sqrt(g * s)).max()
    return out


def orders(errs):
    return [None] + [np.log(errs[i] / errs[i + 1]) / np.log(2) for i in range(len(errs) - 1)]


def discrete_residual(n, scheme, prof, C=-1.0, direction="fwd", L=1.0):
    """Residual of the discrete MMG operator evaluated at the exact manufactured state,
    split into advective, hydrostatic and continuity parts."""
    dx = L / n
    x = np.arange(1, n + 1) * dx
    s = prof.s(x); ue = C / s; b = s + 0.5 * ue**2
    h = s; e = h - b; u = ue
    m, l, r = np.arange(1, n - 1), np.arange(0, n - 2), np.arange(2, n)
    if direction == "fwd":
        adv = (0.25 * ((u[r] + u[m]) * u[r] - (u[m] + u[l]) * u[m]) if scheme == "pa" else u[m] * (u[r] - u[m])) / dx
        Fp, Fm = 0.5 * (u[r] + u[m]) * h[r], 0.5 * (u[m] + u[l]) * h[m]
    else:
        adv = (0.25 * ((u[r] + u[m]) * u[m] - (u[m] + u[l]) * u[l]) if scheme == "pa" else u[m] * (u[m] - u[l])) / dx
        Fp, Fm = 0.5 * (u[r] + u[m]) * h[m], 0.5 * (u[m] + u[l]) * h[l]
    grav = 0.5 * g * (e[r] - e[l]) / dx
    # exact continuous values
    adv_ex = (ue * (-C * prof.ds(x) / s**2))[m]
    grav_ex = -adv_ex
    return dict(Ru=np.abs(adv + grav).max(), Rh=np.abs((Fp - Fm) / dx).max(),
                adv_trunc=np.abs(adv - adv_ex).max(), grav_trunc=np.abs(grav - grav_ex).max())
