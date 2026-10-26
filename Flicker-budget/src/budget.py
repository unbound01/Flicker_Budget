"""
Frequency-noise (FN) budget of a PLL-stabilised qubit LO with a 1/f^3 term.
Baseline = Patra et al., JSSC 53(1), 2018, Eqs. 7-9; extension = flicker (1/f^3) corner.

Model (Patra JSSC 2018, Eqs. 7-9; extension = 1/f^3 term)
    L_osc(df) = alpha/df^2 * (1 + fc/df)          [1/Hz, SSB]
    PLL: brick wall at fPLL, flat in-band floor
        fixed floor      : alpha/fPLL^2
        continuous floor : alpha (1+fc/fPLL)/fPLL^2
    FN^2 = 2 * int_{fa}^{fb} df^2 L(df) d(df)      [Hz^2]

Closed form (both variants, xi = 1-(fa/fPLL)^3):
    FN^2 = 2 alpha [A] + 2 K [B],  K = alpha*fc
    A = (fb-fPLL) + xi*fPLL/3
    B = ln(fb/fPLL) + (xi/3 if continuous else 0)
=> the budget depends on (alpha, K) only; fc_max = K_max(alpha)/alpha.

Run:  python src/budget.py [--no-figures] [--out figures/budget_figures.png]
"""
import argparse
import os

import numpy as np
from scipy.integrate import quad

# ---- Patra parameters ------------------------------------------------------
FA, FB, FPLL = 8.3e3, 10e6, 300e3       # Hz
FN_MAX = 1.9e3                          # Hz rms
ENBW_FACTOR = np.pi**2 / 4              # van Dijk Eq. 10, theta = pi (ENBW = 2.47 fR)


# ---- closed form -----------------------------------------------------------
def _xi(fa, fpll):
    return 1.0 - (fa / fpll) ** 3


def AB(fpll=FPLL, fb=FB, fa=FA, continuous=False):
    if not (fa < fpll < fb):
        raise ValueError(f"PLL bandwidth {fpll:g} Hz must lie between the integration limits {fa:g} and {fb:g} Hz")
    xi = _xi(fa, fpll)
    A = (fb - fpll) + xi * fpll / 3.0
    B = np.log(fb / fpll) + (xi / 3.0 if continuous else 0.0)
    return A, B


def fn2(alpha, fc=0.0, fpll=FPLL, fb=FB, fa=FA, continuous=False):
    A, B = AB(fpll, fb, fa, continuous)
    return 2.0 * alpha * A + 2.0 * alpha * fc * B


def fn(alpha, fc=0.0, **kw):
    return np.sqrt(fn2(alpha, fc, **kw))


def alpha_max(fc=0.0, fn_max=FN_MAX, **kw):
    A, B = AB(**kw)
    return fn_max**2 / (2.0 * (A + fc * B))


def fc_max(alpha, fn_max=FN_MAX, **kw):
    """Largest tolerable corner; <0 means alpha alone already violates the budget."""
    A, B = AB(**kw)
    return (fn_max**2 / (2.0 * alpha) - A) / B


def fn2_free_running(alpha, fc=0.0, fb=FB, fa=FA):
    """No PLL: integral of alpha/df^2*(1+fc/df)*df^2 from fa to fb."""
    return 2.0 * alpha * ((fb - fa) + fc * np.log(fb / fa))


# ---- numerical unit test (NOT a result) -----------------------------------
def fn2_numeric(alpha, fc=0.0, fpll=FPLL, fb=FB, fa=FA, continuous=False):
    floor = alpha * (1 + (fc / fpll if continuous else 0.0)) / fpll**2
    inband = quad(lambda f: f**2 * floor, fa, fpll, epsabs=0, epsrel=1e-12)[0]
    out = quad(lambda f: f**2 * alpha / f**2 * (1 + fc / f), fpll, fb,
               epsabs=0, epsrel=1e-12)[0]
    return 2.0 * (inband + out)


def unit_test():
    worst = 0.0
    for cont in (False, True):
        for alpha in (0.05, 0.184, 0.58):
            for fc in (0.0, 1e5, 5e5, 3e6):
                for fpll in (1e5, 3e5, 1e6):
                    cf = fn2(alpha, fc, fpll=fpll, continuous=cont)
                    nu = fn2_numeric(alpha, fc, fpll=fpll, continuous=cont)
                    worst = max(worst, abs(cf / nu - 1))
    assert worst < 1e-3, worst
    return worst


# ---- phase-noise helpers ---------------------------------------------------
def dbc(x):
    return 10 * np.log10(x)


def ok(cond):
    return "PASS" if cond else "FAIL"


def checkpoints():
    print("=" * 72)
    print("CHECKPOINT 1  baseline (fPLL=300 kHz, fa=8.3 kHz, fb=10 MHz)")
    a = alpha_max()
    print(f"  alpha_max            = {a:.4f} Hz^2/Hz   (expect 0.1842)  {ok(abs(a-0.1842)<5e-4)}")
    inb, at10 = dbc(a / FPLL**2), dbc(a / FB**2)
    print(f"  in-band PN           = {inb:.1f} dBc/Hz     (expect -116.9)  {ok(abs(inb+116.9)<0.1)}")
    print(f"  PN at 10 MHz         = {at10:.1f} dBc/Hz     (expect -147.3)  {ok(abs(at10+147.3)<0.1)}")
    print(f"  Patra alpha=0.2 -> FN = {fn(0.2)/1e3:.3f} kHz (Patra quotes <0.2; exact bound is tighter)")
    print(f"  alpha=0.2 in-band PN = {dbc(0.2/FPLL**2):.1f} dBc/Hz  vs Patra -115 -> 2 dB gap NOT explained by rounding")
    w = unit_test()
    print(f"  closed form vs quad  : worst rel. diff = {w:.1e}  {ok(w<1e-3)}  (unit test, not a result)")

    print("-" * 72)
    print("CHECKPOINT 3  corner tolerance fc_max (kHz)")
    print("  alpha    fixed floor   continuous floor")
    for al in (0.17, 0.15, 0.10):
        print(f"  {al:5.2f}   {fc_max(al)/1e3:9.1f}     {fc_max(al, continuous=True)/1e3:9.1f}")
    print("  expect fixed: 233.2 / 636.9 / 2353   continuous (alpha=0.17): ~213")
    f1, f2, f3 = (fc_max(x) / 1e3 for x in (0.17, 0.15, 0.10))
    print(f"  fixed floor check    {ok(abs(f1-233.2)<1 and abs(f2-636.9)<1 and abs(f3-2353)<5)}")
    fn500 = fn(a, 500e3) / 1e3
    print(f"  FN at alpha_max, fc=500 kHz = {fn500:.3f} kHz (expect 2.063)  {ok(abs(fn500-2.063)<2e-3)}")
    tot = fn2(a, 500e3)
    A, B = AB()
    xi = _xi(FA, FPLL)
    white = 2 * a * (FB - FPLL) / tot
    flick = 2 * a * 500e3 * B / tot
    inband = 2 * a * xi * FPLL / 3 / tot
    print(f"  FN^2 split: white out-of-band {white:.0%}, flicker {flick:.0%}, in-band {inband:.1%}"
          "  (expect ~84/15/0.9)")

    print("-" * 72)
    print("CHECKPOINT 4  f_b sweep at fPLL=300 kHz (headline); fR = fb/ENBW_FACTOR")
    print("  fb[MHz]  alpha_max   equivalent fR[MHz]")
    for fb in (1e6, 2e6, 5e6, 10e6):
        print(f"  {fb/1e6:6.0f}   {alpha_max(fb=fb):8.2f}   {fb/ENBW_FACTOR/1e6:8.2f}")
    print("  expect 2.26 / 1.00 / 0.38 / 0.18")
    got = [alpha_max(fb=f) for f in (1e6, 2e6, 5e6, 10e6)]
    print(f"  {ok(np.allclose(got, [2.26, 1.00, 0.38, 0.18], atol=0.006))}")

    print("-" * 72)
    print("BRIDGE  Gong 2022 measured oscillator (4.2 K): beta=0.041, fC=1.3 MHz")
    for cont in (False, True):
        tag = "continuous" if cont else "fixed     "
        print(f"  {tag}: FN = {fn(0.041, 1.3e6, continuous=cont)/1e3:.3f} kHz, "
              f"fc_max = {fc_max(0.041, continuous=cont)/1e6:.2f} MHz "
              f"(margin {fc_max(0.041, continuous=cont)/1.3e6:.1f}x)")
    print("  NOTE: alpha is at the oscillator carrier (4.65 GHz); x N^2 if multiplied to the qubit LO.")

    print("-" * 72)
    print("GAP 2  Patra: free-running FN = 3.4 kHz (8.3 kHz..10 MHz) vs 'PLL <= 1 MHz easily reaches 1.9 kHz'")
    fnf = 3.4e3
    a_eff = fnf**2 / (2 * (FB - FA))
    print(f"  pure 1/f^2: alpha_eff = {a_eff:.2f} Hz^2/Hz ({a_eff/alpha_max():.1f}x over budget)")
    La = np.log(FB / FA)
    for fpll in (3e5, 1e6):
        for cont in (False, True):
            A, B = AB(fpll=fpll, continuous=cont)
            # fix free-running FN=3.4 kHz: K = (FN^2/2 - alpha (fb-fa))/ln(fb/fa); scan alpha
            al = np.linspace(0, a_eff, 2001)
            K = (fnf**2 / 2 - al * (FB - FA)) / La
            f2 = 2 * al * A + 2 * K * B
            i = np.argmin(f2)
            print(f"  fPLL={fpll/1e3:5.0f} kHz {'cont ' if cont else 'fixed'}: "
                  f"FN range {np.sqrt(f2.min())/1e3:.2f}..{np.sqrt(f2.max())/1e3:.2f} kHz "
                  f"(min at alpha={al[i]:.2f}, i.e. nearly pure 1/f^3)")
    print("  -> for fPLL=1 MHz even the best case sits at ~1.9-2.1 kHz: 'easily' is marginal at best.")
    for fa in (8.3e3, 83e3):
        k = (fnf**2 / 2 - 0.18 * (FB - fa)) / np.log(FB / fa)
        print(f"  fa={fa/1e3:5.1f} kHz -> fc (alpha=0.18) = {k/0.18/1e6:.2f} MHz  (cutoff sensitivity)")
    print("=" * 72)


# ---- figures ---------------------------------------------------------------
def figures(path="figures/budget_figures.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))

    # (a) f_b sweep with f_R axis
    a0 = ax[0, 0]
    fb = np.logspace(6, 7.3, 200)
    a0.loglog(fb / 1e6, [alpha_max(fb=x) for x in fb], "k", label=r"$\alpha_{max}$ (fPLL=300 kHz)")
    a0.axhline(0.041, color="C2", ls="--", label=r"Gong 4.2 K: $\beta$=0.041")
    a0.axhline(3.4e3**2 / (2 * (FB - FA)), color="C3", ls="--", label=r"Patra free-run 3.4 kHz: $\alpha_{eff}$=0.58")
    for x in (1e6, 2e6, 5e6, 10e6):
        a0.plot(x / 1e6, alpha_max(fb=x), "ko")
    a0.set_xlabel("$f_b$ [MHz]")
    a0.set_ylabel(r"$\alpha_{max}$ [Hz$^2$/Hz]")
    a0.set_title("(a) upper band limit is the main lever")
    sec = a0.secondary_xaxis("top", functions=(lambda x: x / ENBW_FACTOR, lambda x: x * ENBW_FACTOR))
    sec.set_xlabel(r"equivalent $f_R$ [MHz] (ENBW = 2.47 $f_R$)")
    a0.legend(fontsize=8)
    a0.grid(True, which="both", alpha=0.3)

    # (b) fc_max vs fPLL
    a1 = ax[0, 1]
    fp = np.logspace(4.7, 6.5, 200)
    for al, c in ((0.17, "C0"), (0.15, "C1"), (0.10, "C2"), (0.041, "C4")):
        y = np.array([fc_max(al, fpll=x) for x in fp])
        y = np.where(y > 0, y, np.nan)
        a1.loglog(fp / 1e3, y / 1e3, c, label=fr"$\alpha$={al}")
        y2 = np.array([fc_max(al, fpll=x, continuous=True) for x in fp])
        a1.loglog(fp / 1e3, np.where(y2 > 0, y2, np.nan) / 1e3, c, ls=":")
    a1.set_xlabel("$f_{PLL}$ [kHz]")
    a1.set_ylabel(r"$f_{c,max}$ [kHz]")
    a1.set_title("(b) corner tolerance (solid: fixed floor, dotted: continuous)")
    a1.legend(fontsize=8)
    a1.grid(True, which="both", alpha=0.3)

    # (c) alpha-fc plane with constant-K hyperbola
    a2 = ax[1, 0]
    al = np.logspace(-2, np.log10(alpha_max() * 0.999), 300)
    a2.loglog(al, fc_max(al) / 1e3, "k", label="fixed floor")
    a2.loglog(al, fc_max(al, continuous=True) / 1e3, "k:", label="continuous floor")
    K = 0.041 * 1.3e6
    a2.loglog(al, K / al / 1e3, "C4--", label=r"const. $K=\alpha f_c$ (Gong 4.2 K)")
    a2.plot(0.041, 1.3e3, "C4*", ms=12, label="Gong: 0.041, 1.3 MHz")
    a2.fill_between(al, fc_max(al) / 1e3, 1e5, color="C3", alpha=0.08, label="budget exceeded (not allowed)")
    a2.set_ylim(10, 1e5)
    a2.set_xlabel(r"$\alpha$ [Hz$^2$/Hz]")
    a2.set_ylabel(r"$f_c$ [kHz]")
    a2.set_title("(c) allowed region (below curve); budget depends on K only")
    a2.legend(fontsize=8)
    a2.grid(True, which="both", alpha=0.3)

    # (d) FN^2 breakdown vs fc at alpha_max-ish
    a3 = ax[1, 1]
    a = alpha_max()
    fcs = np.linspace(0, 1e6, 6)
    A, B = AB()
    xi = _xi(FA, FPLL)
    white = np.full_like(fcs, 2 * a * (FB - FPLL))
    inb = np.full_like(fcs, 2 * a * xi * FPLL / 3)
    fl = 2 * a * fcs * B
    a3.bar(fcs / 1e3, white, width=60, label="white out-of-band")
    a3.bar(fcs / 1e3, inb, width=60, bottom=white, label="in-band")
    a3.bar(fcs / 1e3, fl, width=60, bottom=white + inb, label="flicker")
    a3.axhline(FN_MAX**2, color="r", ls="--", label="budget")
    a3.set_xlabel(r"$f_c$ [kHz] at $\alpha_{max}$")
    a3.set_ylabel(r"FN$^2$ [Hz$^2$]")
    a3.set_title("(d) FN$^2$ breakdown (fixed floor)")
    a3.legend(fontsize=8)

    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=150)
    print(f"saved {path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--out", default="figures/budget_figures.png")
    args = ap.parse_args()
    checkpoints()
    if not args.no_figures:
        figures(args.out)
