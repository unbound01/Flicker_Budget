#!/usr/bin/env python3
"""Phase 2: is the fitted flicker corner identifiable from the digitized 4 K phase-noise curve?

Model:  L(f) = alpha/f^2 * (1 + fc/f)   [fitted in dB, parameters ln(alpha), ln(fc)]

Studies (primary fit range 200 kHz - 10 MHz):
  1. fit over a grid of lower/upper range edges
  2. residuals in dB versus offset
  3. linearised covariance, chi^2 surface and profile likelihood in (alpha, fc)
  4. Monte Carlo propagation of digitization error (independent, common-offset, tilt)
  5. comparison with the independently digitized corner of Patra's Fig. 24
  6. alpha-fc trade-off: correlation, ridge slope, uncertainty of K = alpha*fc

Needs data/patra_digitized.json (generate locally with digitize_patra_figs.py).
Usage:  python scripts/fit_identifiability.py [--data PATH] [--fig figures/fit_sensitivity.png] [--n-mc 2000]
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.optimize import brentq, least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import budget as b  # noqa: E402

PRIMARY = (2e5, 1e7)
SIG_DIG = 0.3          # dB, estimated digitization error per point (docs/validation.md)


def model_db(f, la, lf):
    return 10 * np.log10(np.exp(la) / f**2 * (1 + np.exp(lf) / f))


def select(A, lo, hi):
    m = (A[:, 0] >= lo) & (A[:, 0] <= hi * 1.001)
    return A[m, 0], A[m, 1]


def fit(f, y, p0=(np.log(0.3), np.log(1e6))):
    r = least_squares(lambda p: model_db(f, *p) - y, p0)
    return r


def profile_chi2(f, y, lf_grid, sigma):
    """chi^2 versus ln(fc) with alpha re-optimised analytically (alpha is an additive constant in dB)."""
    out = []
    for lf in lf_grid:
        g = 10 * np.log10(1 + np.exp(lf) / f) - 20 * np.log10(f)
        a_db = np.mean(y - g)                      # 10 log10(alpha*)
        out.append(np.sum((g + a_db - y) ** 2) / sigma**2)
    return np.array(out)


def derived(al, fc):
    fn_free = np.sqrt(b.fn2_free_running(al, fc))
    fn1 = b.fn(al, fc, fpll=1e6)
    fn1c = b.fn(al, fc, fpll=1e6, continuous=True)
    try:
        need = brentq(lambda x: b.fn(al, fc, fpll=x) - 1.9e3, 1e5, 9e6)
    except ValueError:
        need = np.nan
    return fn_free, fn1, fn1c, need


def pct(x, q=(16, 50, 84)):
    return np.nanpercentile(np.asarray(x), q)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "..", "data", "patra_digitized.json"))
    ap.add_argument("--fig", default=os.path.join(HERE, "..", "figures", "fit_sensitivity.png"))
    ap.add_argument("--n-mc", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if not os.path.exists(a.data):
        sys.exit("data file not found; generate it with scripts/digitize_patra_figs.py (see data/README.md)")
    D = json.load(open(a.data))
    A = np.array(D["fig22_PN_dBc_per_Hz_vs_offset_Hz"]["4K"])
    fc24 = D["fig24_corner_Hz_vs_T_K"][0]["corner_Hz"]
    rng = np.random.default_rng(a.seed)

    # ---------------- primary fit ---------------------------------------------------------
    f, y = select(A, *PRIMARY)
    r = fit(f, y)
    la0, lf0 = r.x
    al0, fc0 = np.exp(la0), np.exp(lf0)
    n = len(f)
    rms = np.sqrt(np.mean(r.fun**2))
    sig_eff = np.sqrt(np.sum(r.fun**2) / (n - 2))
    print(f"PRIMARY {PRIMARY[0]/1e3:.0f} kHz-{PRIMARY[1]/1e6:.0f} MHz: n={n} alpha={al0:.4f} fc={fc0/1e6:.3f} MHz "
          f"K={al0*fc0:.3g} Hz^2 rms={rms:.2f} dB  sigma_eff={sig_eff:.2f} dB (vs digitization {SIG_DIG} dB)")
    print(f"  chi2_red with sigma={SIG_DIG} dB: {np.sum(r.fun**2)/SIG_DIG**2/(n-2):.1f}  -> model misfit exceeds digitization error")
    share = (np.log10(fc0) - np.log10(PRIMARY[0])) / (np.log10(PRIMARY[1]) - np.log10(PRIMARY[0]))
    print(f"  flicker term exceeds white term for f < fc: that is {100*min(share,1):.0f}% of the fitted log-frequency span")

    # ---------------- 1+2: range grid ----------------------------------------------------
    lowers, uppers = (1.5e5, 2e5, 3e5, 5e5), (5e6, 1e7, 2e7, 2.9e7)
    grid = {}
    print("\nRANGE GRID: alpha / fc[MHz] / rms[dB] / FN@1MHz(fixed)[kHz]")
    for up in uppers:
        row = []
        for lo in lowers:
            ff, yy = select(A, lo, up)
            rr = fit(ff, yy)
            al, fc = np.exp(rr.x)
            grid[(lo, up)] = (al, fc, np.sqrt(np.mean(rr.fun**2)), b.fn(al, fc, fpll=1e6) / 1e3, len(ff))
            row.append(f"{al:5.3f}/{fc/1e6:4.1f}/{grid[(lo, up)][2]:.2f}/{grid[(lo, up)][3]:.2f}")
        print(f"  upper {up/1e6:4.1f} MHz | " + "  ".join(f"lo {lo/1e3:3.0f}k: {c}" for lo, c in zip(lowers, row)))
    G = np.array([[v[0], v[1], v[3]] for v in grid.values()])
    print(f"  alpha range {G[:,0].min():.3f}-{G[:,0].max():.3f}; fc range {G[:,1].min()/1e6:.2f}-{G[:,1].max()/1e6:.2f} MHz; "
          f"FN@1MHz {G[:,2].min():.2f}-{G[:,2].max():.2f} kHz")
    Kg = G[:, 0] * G[:, 1]
    print(f"  K=alpha*fc range {Kg.min():.3g}-{Kg.max():.3g} Hz^2 (factor {Kg.max()/Kg.min():.2f}); alpha spans a factor "
          f"{G[:,0].max()/G[:,0].min():.1f}, fc a factor {G[:,1].max()/G[:,1].min():.1f}")

    # ---------------- 3: linearised covariance, chi2 surface, profile ----------------------
    J = r.jac
    cov_unit = np.linalg.inv(J.T @ J)
    print("\nLINEARISED UNCERTAINTY (log parameters)")
    for name, s in (("digitization-only sigma=0.3 dB", SIG_DIG), ("inflated sigma=residual rms", sig_eff)):
        cov = cov_unit * s**2
        sa, sf = np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1])
        corr = cov[0, 1] / (sa * sf)
        sK = np.sqrt(cov[0, 0] + cov[1, 1] + 2 * cov[0, 1])
        print(f"  {name}: alpha x/÷ {np.exp(sa):.2f}, fc x/÷ {np.exp(sf):.2f}, K x/÷ {np.exp(sK):.2f}; corr(ln a, ln fc)={corr:+.3f}")
    cov = cov_unit * sig_eff**2
    slope = cov[0, 1] / cov[0, 0]
    print(f"  ridge slope d ln fc / d ln alpha = {slope:+.2f} (constant K would be -1.00); cond(J)={np.linalg.cond(J):.1f}")

    lf_grid = np.log(np.logspace(np.log10(8e5), np.log10(3e7), 400))
    prof = profile_chi2(f, y, lf_grid, sig_eff)
    prof -= prof.min()
    fcs = np.exp(lf_grid)
    def interval(th):
        ok = fcs[prof <= th]
        return ok.min(), ok.max()
    lo68, hi68 = interval(1.0)
    lo95, hi95 = interval(3.84)
    print(f"\nPROFILE LIKELIHOOD (sigma_eff): fc 68% [{lo68/1e6:.2f}, {hi68/1e6:.2f}] MHz, 95% [{lo95/1e6:.2f}, {hi95/1e6:.2f}] MHz")
    d24 = np.interp(np.log(fc24), lf_grid, prof)
    print(f"  Fig. 24 corner {fc24/1e6:.2f} MHz: delta chi2 = {d24:.2f} relative to best fit "
          f"({'inside' if d24 <= 3.84 else 'outside'} the 95% interval, {'inside' if d24 <= 1 else 'outside'} 68%)")
    g = 10 * np.log10(1 + fc24 / f) - 20 * np.log10(f)
    a24 = 10 ** (np.mean(y - g) / 10)
    print(f"  refit with fc fixed to Fig. 24 value: alpha={a24:.3f}, rms={np.sqrt(np.mean((g+np.mean(y-g)-y)**2)):.2f} dB, "
          f"K={a24*fc24:.3g} Hz^2 (best fit K={al0*fc0:.3g})")

    # chi2 surface for the figure
    la_ax = np.linspace(np.log(0.03), np.log(0.45), 160)
    lf_ax = np.linspace(np.log(6e5), np.log(3e7), 160)
    LA, LF = np.meshgrid(la_ax, lf_ax)
    res = model_db(f[:, None, None], LA[None], LF[None]) - y[:, None, None]
    chi2 = np.sum(res**2, axis=0) / sig_eff**2
    chi2 -= chi2.min()

    # ---------------- 4: Monte Carlo of digitization error ---------------------------------
    scen = {
        "A independent 0.3 dB": lambda: rng.normal(0, 0.3, n),
        "B common offset 0.3 dB": lambda: np.full(n, rng.normal(0, 0.3)),
        "C tilt 1 dB at 29 MHz": lambda: rng.normal(0, 1.0) * (np.log10(f) - 6) / (np.log10(2.9e7) - 6),
        "D all of A+B+C": lambda: (rng.normal(0, 0.2, n) + rng.normal(0, 0.2)
                                   + rng.normal(0, 1.0) * (np.log10(f) - 6) / (np.log10(2.9e7) - 6)),
    }
    mc = {}
    print(f"\nMONTE CARLO ({a.n_mc} refits per scenario): 16 / 50 / 84 percentiles")
    for name, gen in scen.items():
        rows = []
        for _ in range(a.n_mc):
            rr = fit(f, y + gen(), p0=(la0, lf0))
            al, fc = np.exp(rr.x)
            rows.append((al, fc, *derived(al, fc)))
        rows = np.array(rows)
        mc[name] = rows
        P = [pct(rows[:, i]) for i in range(6)]
        PK = pct(rows[:, 0] * rows[:, 1])
        print(f"  {name:24s} K[Hz^2] {PK[0]:.3g}/{PK[1]:.3g}/{PK[2]:.3g}")
        print(f"  {name:24s} alpha {P[0][0]:.3f}/{P[0][1]:.3f}/{P[0][2]:.3f} | fc[MHz] {P[1][0]/1e6:.2f}/{P[1][1]/1e6:.2f}/{P[1][2]/1e6:.2f} | "
              f"FNfree[kHz] {P[2][0]/1e3:.2f}/{P[2][1]/1e3:.2f}/{P[2][2]/1e3:.2f} | FN@1M fixed {P[3][0]/1e3:.2f}/{P[3][1]/1e3:.2f}/{P[3][2]/1e3:.2f} "
              f"| PLL BW for 1.9 kHz [MHz] {P[5][0]/1e6:.2f}/{P[5][1]/1e6:.2f}/{P[5][2]/1e6:.2f} | P(FN@1M<1.9)={np.mean(rows[:,3]<1.9e3):.3f}")

    # ---------------- figure ---------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    # (a) residuals vs offset
    for (lo, up), c in zip(((1.5e5, 2.9e7), (2e5, 1e7), (3e5, 1e7), (5e5, 1e7)), ("C0", "C3", "C2", "C1")):
        ff, yy = select(A, lo, up)
        rr = fit(ff, yy)
        ax[0, 0].semilogx(A[:, 0], model_db(A[:, 0], *rr.x) - A[:, 1], c + "o-", ms=3, lw=1,
                          label=f"fit {lo/1e3:.0f} kHz-{up/1e6:.1f} MHz")
    ax[0, 0].axhspan(-SIG_DIG, SIG_DIG, color="gray", alpha=0.2, label="±0.3 dB digitization")
    ax[0, 0].axhline(0, color="k", lw=0.5)
    ax[0, 0].set_xlabel("offset [Hz]"); ax[0, 0].set_ylabel("model − data [dB]")
    ax[0, 0].set_title("(a) residuals vs offset (all fits, full range)"); ax[0, 0].legend(fontsize=7)
    # (b) chi2 contours
    cs = ax[0, 1].contour(np.exp(LA), np.exp(LF) / 1e6, chi2, levels=[2.30, 6.18, 11.8], colors=["C3", "C1", "C0"])
    ax[0, 1].clabel(cs, fmt={2.30: "68%", 6.18: "95%", 11.8: "99.7%"}, fontsize=7)
    ax[0, 1].plot(al0, fc0 / 1e6, "k*", ms=12, label="primary fit")
    ax[0, 1].plot(G[:, 0], G[:, 1] / 1e6, "k.", ms=5, label="range-grid fits")
    aa = np.logspace(np.log10(0.03), np.log10(0.45), 50)
    ax[0, 1].plot(aa, al0 * fc0 / aa / 1e6, "g--", lw=1, label="constant K = α·f_c")
    ax[0, 1].axhline(fc24 / 1e6, color="m", ls=":", label=f"Fig. 24: {fc24/1e6:.1f} MHz")
    ax[0, 1].set_xscale("log"); ax[0, 1].set_yscale("log")
    ax[0, 1].set_xlabel("α [Hz²/Hz]"); ax[0, 1].set_ylabel("f_c [MHz]")
    ax[0, 1].set_title("(b) Δχ² surface (σ = residual rms)"); ax[0, 1].legend(fontsize=7)
    # (c) profile
    ax[0, 2].semilogx(fcs / 1e6, prof, "k")
    ax[0, 2].axhline(1, color="C3", ls="--", lw=0.8); ax[0, 2].axhline(3.84, color="C1", ls="--", lw=0.8)
    ax[0, 2].axvline(fc0 / 1e6, color="k", ls=":"); ax[0, 2].axvline(fc24 / 1e6, color="m", ls=":", label="Fig. 24 corner")
    ax[0, 2].set_ylim(0, 12); ax[0, 2].set_xlabel("f_c [MHz]"); ax[0, 2].set_ylabel("Δχ² (α re-optimised)")
    ax[0, 2].set_title("(c) profile likelihood in f_c"); ax[0, 2].legend(fontsize=7)
    # (d) range grid
    for up in uppers:
        ax[1, 0].plot([l / 1e3 for l in lowers], [grid[(l, up)][1] / 1e6 for l in lowers], "o-", label=f"upper {up/1e6:.1f} MHz")
    ax[1, 0].axhline(fc24 / 1e6, color="m", ls=":", label="Fig. 24")
    ax[1, 0].set_xlabel("lower fit edge [kHz]"); ax[1, 0].set_ylabel("fitted f_c [MHz]")
    ax[1, 0].set_title("(d) fitted corner vs fit range"); ax[1, 0].legend(fontsize=7)
    # (e) MC scatter
    for (name, rows), c in zip(mc.items(), ("C0", "C1", "C2", "C3")):
        if name.startswith(("A", "D")):
            ax[1, 1].plot(rows[::5, 0], rows[::5, 1] / 1e6, ".", ms=3, alpha=0.5, color=c, label=name)
    ax[1, 1].plot(aa, al0 * fc0 / aa / 1e6, "g--", lw=1, label="constant K")
    ax[1, 1].set_xscale("log"); ax[1, 1].set_yscale("log"); ax[1, 1].set_ylim(1.5, 40)
    ax[1, 1].set_xlabel("α [Hz²/Hz]"); ax[1, 1].set_ylabel("f_c [MHz]")
    ax[1, 1].set_title("(e) digitization-error Monte Carlo"); ax[1, 1].legend(fontsize=7)
    # (f) FN@1MHz
    for (name, rows), c in zip(mc.items(), ("C0", "C1", "C2", "C3")):
        ax[1, 2].hist(rows[:, 3] / 1e3, bins=40, histtype="step", color=c, label=name)
    ax[1, 2].axvline(1.9, color="k", ls="--", label="1.9 kHz budget")
    ax[1, 2].set_xlabel("FN with 1 MHz PLL, fixed floor [kHz]"); ax[1, 2].set_ylabel("samples")
    ax[1, 2].set_title("(f) predicted FN, 1 MHz PLL"); ax[1, 2].legend(fontsize=7)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(a.fig)), exist_ok=True)
    fig.savefig(a.fig, dpi=140)
    print("\nsaved", a.fig)


if __name__ == "__main__":
    main()
