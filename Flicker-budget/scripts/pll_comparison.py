#!/usr/bin/env python3
"""Integrated RMS frequency noise versus PLL bandwidth for three models, identical limits (8.3 kHz - 10 MHz).

  Case 1  Patra's original approximation: oscillator = alpha/f^2 only (alpha from the measured PN at 10 MHz),
          brick-wall PLL, flat floor equal to the oscillator at f_PLL.
  Case 2  flicker-aware brick wall: fitted alpha/f^2 (1 + fc/f); floor "fixed" = alpha/f_PLL^2 (Patra's formula)
          or "continuous" = alpha (1 + fc/f_PLL)/f_PLL^2 (full oscillator PN at f_PLL).
  Case 3  generic linear phase-domain PLL (src/pll.py):  L_out = |H_lp|^2 L_in + |1-H_lp|^2 L_osc,
          generic first-order, type-II and second-order low-pass loops (NOT circuit-specific models), plus
          'Gong-match' = type-II with damping 1.6 (Gong's loop type and damping, without its C_P pole or ZOH).
          In-band floor scenarios, kept separate:
            white      L_in = alpha/f3^2          Patra's formula       (like-for-like with Case 2 "fixed")
            flicker    L_in = L_osc(f3)           full spectrum at f3   (like-for-like with Case 2 "continuous")
            absolute   flat levels in dBc/Hz at 6.3 GHz: -115 (Patra's quoted number) and the SCENARIOS -122, -135, -140
            zero       noiseless input (thought experiment)

The -122 dBc/Hz scenario scales Gong 2023's -130 dBc/Hz at the 2.5 GHz divided output by 20*log10(6.3/2.5);
it is an illustrative scenario, not a demonstrated input-referred floor of the modelled PLL.
No carrier scaling is applied to the oscillator: the 6.3 GHz spectrum is treated as the target oscillator.

Needs the locally generated data/patra_digitized.json (see data/README.md).
Usage: python scripts/pll_comparison.py [--data PATH] [--fig1 PATH] [--fig2 PATH]
"""
import argparse
import json
import os
import sys

import numpy as np
import scipy
from scipy.optimize import brentq, least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import budget as b  # noqa: E402
import pll  # noqa: E402

BUDGET = 1.9e3
SCALED_FLOOR_DB = -130 + 20 * np.log10(6.3 / 2.5)   # = -122.0: scenario, see docstring
PATRA_FLOOR_DB = -115.0                             # Patra's quoted in-band requirement
BW = np.logspace(5, 6.95, 70)                       # closed-loop -3 dB bandwidth (or brick-wall transition)
LOOPS = {"first-order": ("first-order", None), "type-II": ("type-II (PI zero)", 1.0),
         "LP2": ("2nd-order low-pass", 1.6),        # generic second-order low-pass, damping 1.6 (not Gong's loop type)
         "Gong-match": pll.GONG_MATCHED}            # type-II, damping 1.6: the generic family closest to Gong's DAPLL
_trapz = getattr(np, "trapezoid", None) or np.trapz


def fit_primary(A, lo=2e5, hi=1e7):
    m = (A[:, 0] >= lo) & (A[:, 0] <= hi * 1.001)
    f, y = A[m, 0], A[m, 1]
    r = least_squares(lambda p: 10 * np.log10(np.exp(p[0]) / f**2 * (1 + np.exp(p[1]) / f)) - y,
                      np.log((0.3, 1e6)))
    return np.exp(r.x)


def floor_value(floor, f3, al, fc, osc=None):
    if floor == "white":
        return al / f3**2
    if floor == "flicker":
        return float(osc(np.array([f3]))[0]) if osc is not None else pll.l_osc(f3, al, fc)
    if floor == "zero":
        return 0.0
    return pll.dbc_to_lin(floor)                    # absolute level in dBc/Hz


def curve(loop, floor, al, fc, osc=None):
    arch, zeta = LOOPS[loop] if isinstance(loop, str) else loop
    return np.array([pll.fn_rms(pll.make_h_lp(arch, f3, zeta), al, fc, floor_value(floor, f3, al, fc, osc), osc=osc)
                     for f3 in BW])


def crossing(x, y, th=BUDGET):
    """Lowest bandwidth at which y falls to th (log interpolation); nan if never."""
    for i in range(len(x) - 1):
        if y[i] > th >= y[i + 1]:
            t = (y[i] - th) / (y[i] - y[i + 1])
            return float(np.exp(np.log(x[i]) + t * (np.log(x[i + 1]) - np.log(x[i]))))
    return np.nan if y.min() > th else float(x[0])


def required_floor(loop, al, fc):
    """Largest flat in-band floor [dBc/Hz] meeting 1.9 kHz at each bandwidth (nan: oscillator alone fails)."""
    arch, zeta = LOOPS[loop]
    res = []
    for f3 in BW:
        h = pll.make_h_lp(arch, f3, zeta)
        g = lambda ll: pll.fn_rms(h, al, fc, 10 ** ll) - BUDGET
        res.append(np.nan if g(-30) > 0 else (-80.0 if g(-8) < 0 else 10 * brentq(g, -30, -8)))
    return np.array(res)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "..", "data", "patra_digitized.json"))
    ap.add_argument("--fig1", default=os.path.join(HERE, "..", "figures", "pll_three_cases.png"))
    ap.add_argument("--fig2", default=os.path.join(HERE, "..", "figures", "pll_sensitivity.png"))
    a = ap.parse_args()
    if not os.path.exists(a.data):
        sys.exit("data file not found; generate it with scripts/digitize_patra_figs.py (see data/README.md)")
    print(f"numpy {np.__version__}, scipy {scipy.__version__}")
    D = json.load(open(a.data))
    A = np.array(D["fig22_PN_dBc_per_Hz_vs_offset_Hz"]["4K"])
    al, fc = fit_primary(A)
    l10 = np.interp(np.log10(1e7), np.log10(A[:, 0]), A[:, 1])
    al_p = 10 ** (l10 / 10) * 1e14
    fr = np.sqrt(b.fn2_free_running(al, fc))
    print(f"oscillator (4 K, primary fit): alpha={al:.4f} Hz^2/Hz, fc={fc/1e6:.2f} MHz; free-running FN={fr/1e3:.2f} kHz")
    print(f"Case 1 alpha from measured PN at 10 MHz ({l10:.2f} dBc/Hz): {al_p:.4f}; scaled-floor scenario {SCALED_FLOOR_DB:.1f} dBc/Hz")

    c1 = np.array([b.fn(al_p, 0.0, fpll=f) for f in BW])
    c2f = np.array([b.fn(al, fc, fpll=f) for f in BW])
    c2c = np.array([b.fn(al, fc, fpll=f, continuous=True) for f in BW])
    floors = ["white", "flicker", SCALED_FLOOR_DB, -135.0, -140.0, PATRA_FLOOR_DB, "zero"]
    c3 = {(lp, fl): curve(lp, fl, al, fc) for lp in LOOPS for fl in floors}

    def name(fl):
        if isinstance(fl, str):
            return {"white": "white-FM floor alpha/f3^2", "flicker": "flicker-aware floor L_osc(f3)", "zero": "perfect (zero) input"}[fl]
        return f"scenario floor {fl:.0f} dBc/Hz"

    def row(label, y):
        i = int(np.argmin(y))
        cr = crossing(BW, y)
        at = [np.interp(np.log(x), np.log(BW), y) / 1e3 for x in (3e5, 1e6, 3e6)]
        print(f"  {label:50s} FN@0.3/1/3 MHz = {at[0]:5.2f}/{at[1]:5.2f}/{at[2]:5.2f} kHz | min {y[i]/1e3:5.2f} kHz at {BW[i]/1e6:4.2f} MHz | "
              f"1.9 kHz at: {('%.2f MHz' % (cr/1e6)) if np.isfinite(cr) else 'never'}")

    print("\nINTEGRATED FN versus bandwidth (8.3 kHz-10 MHz), fitted oscillator spectrum:")
    row("Case 1 Patra (alpha only)", c1)
    row("Case 2 brick wall, fixed floor (alpha/f^2)", c2f)
    row("Case 2 brick wall, continuous floor (L_osc)", c2c)
    bw_zero = np.array([np.sqrt(2 * al * ((b.FB - f) + fc * np.log(b.FB / f))) for f in BW])
    row("Case 2 brick wall, perfect input", bw_zero)
    for lp in LOOPS:
        for fl in floors:
            row(f"Case 3 {lp}, {name(fl)}", c3[(lp, fl)])

    # ---- direct digitized spectrum instead of the two-term fit --------------------------------------
    osc_d = pll.l_osc_from_data(A[:, 0], A[:, 1], al, fc)
    f_g = np.logspace(np.log10(b.FA), np.log10(b.FB), 200001)
    fr_d = np.sqrt(2 * _trapz(osc_d(f_g) * f_g**3, np.log(f_g)))
    print(f"\nDIRECT DIGITIZED SPECTRUM (log-log interpolation above {A[0,0]/1e3:.0f} kHz, fit below): free-running FN {fr_d/1e3:.3f} kHz "
          f"(two-term fit {fr/1e3:.3f} kHz)")
    print("  loop / floor -> FN at 1 MHz and bandwidth for 1.9 kHz")
    direct = {}
    for lp in LOOPS:
        parts = []
        for fl in ("white", "flicker", SCALED_FLOOR_DB, "zero"):
            y = curve(lp, fl, al, fc, osc=osc_d)
            direct[(lp, fl)] = y
            cr = crossing(BW, y)
            parts.append(f"{np.interp(np.log(1e6), np.log(BW), y)/1e3:.2f} kHz / {('%.2f MHz' % (cr/1e6)) if np.isfinite(cr) else 'never'}")
        print(f"  {lp:12s} white {parts[0]} | measured-PN floor {parts[1]} | {SCALED_FLOOR_DB:.0f} dBc/Hz {parts[2]} | zero {parts[3]}")

    print("\nMATCHED FLOOR levels at the loop bandwidth [dBc/Hz]: white-FM alpha/f3^2 | flicker-aware L_osc(f3)")
    for f3 in (3e5, 1e6, 3e6):
        print(f"  {f3/1e6:.1f} MHz: {10*np.log10(al/f3**2):.1f} | {10*np.log10(pll.l_osc(f3, al, fc)):.1f}")

    # ---- oscillator-fit sensitivity: crossing computed for each of the 16 fit windows ------------------
    fits = [fit_primary(A, lo, hi) for lo in (1.5e5, 2e5, 3e5, 5e5) for hi in (5e6, 1e7, 2e7, 2.9e7)]
    env, cross = {}, {}
    for lp, fl in (("LP2", "white"), ("LP2", "flicker"), ("type-II", "white"), ("first-order", "white"), ("LP2", SCALED_FLOOR_DB)):
        cs = [curve(lp, fl, a_i, f_i) for a_i, f_i in fits]
        env[(lp, fl)] = (np.min(cs, axis=0), np.max(cs, axis=0))
        cross[(lp, fl)] = [crossing(BW, c) for c in cs]
    print("\nOSCILLATOR-FIT SENSITIVITY: bandwidth for 1.9 kHz computed per fit window (16 windows; a sensitivity spread, not a confidence interval)")
    for k, v in cross.items():
        v = np.array(v)
        ok = v[np.isfinite(v)]
        print(f"  {k[0]:12s} {name(k[1]):32s}: " + (f"{ok.min()/1e6:.2f}-{ok.max()/1e6:.2f} MHz ({len(ok)}/16 windows reach it)" if len(ok) else "never reached in any window"))

    rf = {lp: required_floor(lp, al, fc) for lp in ("type-II", "LP2")}
    print("\nREQUIRED in-band floor to meet 1.9 kHz [dBc/Hz] at given bandwidth (nan: oscillator alone fails for that loop):")
    for k, y in rf.items():
        vals = [np.interp(np.log(x0), np.log(BW), y) for x0 in (1e6, 3e6, 5e6)]
        j = int(np.nanargmax(y))
        print(f"  {k:8s}: " + ", ".join(f"{v:.1f} @ {x0/1e6:g} MHz" for v, x0 in zip(vals, (1e6, 3e6, 5e6))) + f"; best {y[j]:.1f} at {BW[j]/1e6:.2f} MHz")

    # ---------------- figure 1: three cases ------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = BW / 1e6
    fig, ax = plt.subplots(figsize=(10, 8.8))
    ax.semilogx(x, c1 / 1e3, "k--", lw=2, label="Case 1: Patra's approximation (α/f² only)")
    ax.semilogx(x, c2f / 1e3, "C0-", lw=2, label="Case 2: flicker-aware brick wall, fixed floor (α/f²)")
    ax.semilogx(x, c2c / 1e3, "C0:", lw=2, label="Case 2: flicker-aware brick wall, continuous floor (L_osc)")
    ax.semilogx(x, c3[("type-II", "white")] / 1e3, "C1-", lw=2, label="Case 3: generic type-II loop, white-FM floor α/f²")
    ax.semilogx(x, c3[("LP2", "white")] / 1e3, "C2-", lw=2, label="Case 3: generic 2nd-order low-pass loop, white-FM floor α/f²")
    ax.semilogx(x, c3[("type-II", "flicker")] / 1e3, "C1-.", lw=2, label="Case 3: generic type-II loop, flicker-aware floor L_osc(f₃)")
    ax.semilogx(x, c3[("LP2", "flicker")] / 1e3, "C2-.", lw=2, label="Case 3: generic 2nd-order low-pass loop, flicker-aware floor L_osc(f₃)")
    ax.semilogx(x, c3[("LP2", SCALED_FLOOR_DB)] / 1e3, "C4--", lw=2,
                label=f"Case 3: 2nd-order low-pass loop, scenario floor −{abs(SCALED_FLOOR_DB):.0f} dBc/Hz (scaled from Gong)")
    ax.semilogx(x, c3[("type-II", "zero")] / 1e3, "C1:", lw=2, label="Case 3: type-II loop, perfect (zero) input noise (thought experiment)")
    ax.axhline(1.9, color="r", lw=1.5, label="1.9 kHz budget")
    ax.axhline(fr / 1e3, color="gray", ls="-.", lw=1, label="free-running oscillator")
    ax.set_ylim(1, 6); ax.set_xlim(x[0], x[-1])
    ax.set_xlabel("PLL bandwidth [MHz] (brick-wall transition, or closed-loop −3 dB)")
    ax.set_ylabel("integrated rms frequency noise, 8.3 kHz–10 MHz [kHz]")
    ax.set_title("Integrated FN versus PLL bandwidth: three models, identical integration limits")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    fig.tight_layout(); os.makedirs(os.path.dirname(os.path.abspath(a.fig1)), exist_ok=True); fig.savefig(a.fig1, dpi=150)

    # ---------------- figure 2: sensitivity -------------------------------------------------------------
    fig, ax = plt.subplots(2, 2, figsize=(14, 9.5))
    for lab, fl, st in (("no floor (perfect input)", "zero", "k:"), ("white-FM floor α/f₃²", "white", "k-"),
                        ("flicker-aware floor L_osc(f₃)", "flicker", "k-."), ("scenario floor −140 dBc/Hz", -140.0, "C0-"),
                        ("scenario floor −135 dBc/Hz", -135.0, "C5-"), ("scenario floor −130 dBc/Hz", -130.0, "C1-"),
                        (f"scenario floor −{abs(SCALED_FLOOR_DB):.0f} dBc/Hz (scaled from Gong)", SCALED_FLOOR_DB, "C2-"),
                        ("floor −115 dBc/Hz (Patra's quoted number)", PATRA_FLOOR_DB, "C3-")):
        ax[0, 0].semilogx(x, curve("type-II", fl, al, fc) / 1e3, st, label=lab)
    ax[0, 0].axhline(1.9, color="r"); ax[0, 0].set_ylim(1, 8)
    ax[0, 0].set_title("(a) in-band floor scenarios, generic type-II loop"); ax[0, 0].legend(fontsize=7)
    for arch, zs, cs in (("type-II (PI zero)", (0.7, 1.0, 1.6), "C1"), ("2nd-order low-pass", (0.7, 1.0, 1.6), "C2")):
        for z, ls in zip(zs, (":", "-", "--")):
            ax[0, 1].semilogx(x, curve((arch, z), "white", al, fc) / 1e3, cs + ls, label=f"{arch}, ζ={z}")
    ax[0, 1].semilogx(x, c3[("first-order", "white")] / 1e3, "k-", label="first-order")
    ax[0, 1].semilogx(x, curve(("4th-order Butterworth (steep)", 1.0), "white", al, fc) / 1e3, "m-", label="4th-order Butterworth (steep)")
    ax[0, 1].axhline(1.9, color="r"); ax[0, 1].set_ylim(1, 8)
    ax[0, 1].set_title("(b) loop family and damping, white-FM floor α/f₃²"); ax[0, 1].legend(fontsize=7)
    for k, c in ((("LP2", "white"), "C2"), (("LP2", SCALED_FLOOR_DB), "C4")):
        lo_c, hi_c = env[k]
        ax[1, 0].fill_between(x, lo_c / 1e3, hi_c / 1e3, color=c, alpha=0.3, label=f"LP2, {name(k[1])}: spread over 16 fit windows")
        ax[1, 0].semilogx(x, c3[k] / 1e3, c + "-")
    ax[1, 0].semilogx(x, direct[("LP2", "white")] / 1e3, "k--", lw=1.5, label="LP2, α/f₃² floor, direct digitized spectrum")
    ax[1, 0].axhline(1.9, color="r"); ax[1, 0].set_ylim(1, 8)
    ax[1, 0].set_title("(c) oscillator sensitivity (spread of fit windows, not a confidence interval)"); ax[1, 0].legend(fontsize=7)
    for k, c in (("type-II", "C1"), ("LP2", "C2")):
        ax[1, 1].semilogx(x, rf[k], c + "-", label=f"max in-band floor, generic {k} loop")
    ax[1, 1].axhline(SCALED_FLOOR_DB, color="g", ls="--", label=f"scenario floor −{abs(SCALED_FLOOR_DB):.0f} (scaled from Gong)")
    ax[1, 1].axhline(PATRA_FLOOR_DB, color="m", ls=":", label="Patra's quoted −115")
    ax[1, 1].set_title("(d) in-band floor required to meet 1.9 kHz (gaps: not achievable for that loop)")
    ax[1, 1].legend(fontsize=7)
    for axx in ax.flat:
        axx.set_xlabel("PLL bandwidth [MHz]"); axx.grid(True, which="both", alpha=0.3)
    for axx in (ax[0, 0], ax[0, 1], ax[1, 0]):
        axx.set_ylabel("integrated FN [kHz]")
    ax[1, 1].set_ylabel("required floor [dBc/Hz, at 6.3 GHz]")
    fig.tight_layout(); fig.savefig(a.fig2, dpi=140)
    print("\nsaved", a.fig1, "and", a.fig2)


if __name__ == "__main__":
    main()
