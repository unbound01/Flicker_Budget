#!/usr/bin/env python3
"""A0: what does Gong's MEASURED 4.2 K PLL output (JSSC 2023, Fig. 13(b)) give for
   (1) the instrument's own integrated numbers (validation of the digitization),
   (2) Gong's own fidelity numbers via his Eqs. 1-2 (validation of the fidelity pipeline),
   (3) Patra's 1.9 kHz rms FN budget over 8.3 kHz - 10 MHz,
   (4) gate fidelity including the measured flicker region (Gong footnote 2)?

Conventions (checked against the instrument box): L is SSB; int noise = int L df; PM = sqrt(2 int L df);
FM = sqrt(2 int f^2 L df); jitter = PM / (2 pi f0). Dividing a carrier by N lowers L by 20log10(N) and FM by N.
Usage: python scripts/a0_gong_measured.py [--data data/gong_fig13b.json] [--data300 data/gong_fig13a.json]
See docs/a0_gong_measured.md.
"""
import argparse
import json
import os

import numpy as np

F_TEST = 2.500018e9          # Hz, analyzer signal frequency (after /4)
N_DIV = 4
F_VCO = N_DIV * F_TEST       # ~10 GHz


def spectrum(data, f_lo_ext=8.0e3, f_hi_ext=3.0e7):
    t = np.array(data["trace"])
    f, d = t[:, 0], t[:, 1]
    lf = np.log10(f)
    # end extensions with the local log-log slope (only 10.0-10.5 kHz and 28.5-30 MHz are outside the trace;
    # 8.3-10 kHz is needed for the budget band)
    sl_lo = np.polyfit(lf[:6], d[:6], 1)
    sl_hi = np.polyfit(lf[-6:], d[-6:], 1)

    def L_db(x):
        x = np.asarray(x, float)
        lx = np.log10(x)
        out = np.interp(lx, lf, d)
        out = np.where(x < f[0], np.polyval(sl_lo, lx), out)
        out = np.where(x > f[-1], np.polyval(sl_hi, lx), out)
        return out
    return L_db, (sl_lo[0], sl_hi[0])


def integ(fun, a, b, n=200001):
    x = np.logspace(np.log10(a), np.log10(b), n)
    y = fun(x) * x                                  # d f = f d ln f
    return np.trapezoid(y, np.log(x))


def h_lo2(f, fr):
    """Gong Eq. 2 (theta = pi), worst-case |H_LO|^2; removable singularity at f = fR handled by a tiny offset."""
    u = np.asarray(f, float) / fr
    u = np.where(np.abs(u - 1) < 1e-9, 1 + 1e-9, u)
    return (1 + np.cos(np.pi * u)) * (1 + u**2) / (2 * (1 - u**2) ** 2)


def infidelity(Llin, fr, a, b):
    return integ(lambda x: 2 * Llin(x) * x**2 / fr**2 * h_lo2(x, fr), a, b, n=400001)


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--data", default=os.path.join(here, "..", "data", "gong_fig13b.json"))
    ap.add_argument("--data300", default=os.path.join(here, "..", "data", "gong_fig13a.json"))
    a = ap.parse_args()
    D = json.load(open(a.data))
    L_db, slopes = spectrum(D)
    Lt = lambda x: 10 ** (L_db(x) / 10)                       # at the 2.5 GHz test output
    Lv = lambda x: 10 ** ((L_db(x) + 20 * np.log10(N_DIV)) / 10)   # referred to the ~10 GHz VCO
    ref = D["instrument_reference"]["integrated_10k_30M"]
    print(f"end slopes used for extension: {slopes[0]:.1f} dB/dec (low), {slopes[1]:.1f} dB/dec (high)")

    print("\n(1) VALIDATION vs instrument box, 10 kHz - 30 MHz, at 2.5 GHz")
    I = integ(Lt, 1e4, 3e7)
    pm = np.sqrt(2 * I)
    fm = np.sqrt(2 * integ(lambda x: Lt(x) * x**2, 1e4, 3e7))
    jit = pm / (2 * np.pi * F_TEST)
    for name, got, want, unit in (("int noise", 10 * np.log10(I), ref["int_noise_dBc"], "dBc"),
                                  ("PM", pm * 1e6, ref["PM_urad"], "urad"),
                                  ("FM", fm, ref["FM_Hz"], "Hz"),
                                  ("jitter", jit * 1e15, ref["jitter_fs"], "fs")):
        dev = got - want if unit == "dBc" else 100 * (got / want - 1)
        print(f"  {name:9s} digitized {got:10.3f} {unit:4s} | instrument {want:10.3f} | "
              + (f"diff {dev:+.2f} dB" if unit == "dBc" else f"diff {dev:+.1f}%"))

    print("\n(2) VALIDATION vs Gong's stated fidelity (Eqs. 1-2, theta = pi, L at the 10 GHz PLL output)")
    for fr, stated in ((1e6, 99.9994), (10e6, 99.9998)):
        for lo, hi in ((1e4, 3e7), (1e3, 1e8)):
            inf = infidelity(Lv, fr, lo, hi)
            print(f"  fR = {fr/1e6:4.0f} MHz, band {lo/1e3:g} kHz-{hi/1e6:g} MHz: F = {100*(1-inf):.5f}%  (1-F = {inf:.2e}); "
                  f"Gong: {stated}% (1-F = {1-stated/100:.1e})")

    print("\n(3) PATRA BUDGET: integrated rms FN over 8.3 kHz - 10 MHz (target 1.9 kHz)")
    fn_t = np.sqrt(2 * integ(lambda x: Lt(x) * x**2, 8.3e3, 1e7))
    rows = (("2.5 GHz test output (as measured)", fn_t),
            ("10 GHz PLL output (x4)", fn_t * N_DIV),
            ("6 GHz, same time noise (x 6/2.5)", fn_t * 6e9 / F_TEST))
    for lab, v in rows:
        print(f"  {lab:38s}: {v/1e3:6.3f} kHz rms  ({'meets' if v <= 1.9e3 else 'exceeds'} 1.9 kHz, ratio {v/1.9e3:.2f})")
    # where does FN^2 come from (10 GHz plane)?
    edges = [8.3e3, 1e5, 1e6, 4e6, 1e7]
    tot = 2 * integ(lambda x: Lv(x) * x**2, 8.3e3, 1e7)
    parts = [2 * integ(lambda x: Lv(x) * x**2, lo, hi) / tot for lo, hi in zip(edges[:-1], edges[1:])]
    print("  FN^2 share by offset band (10 GHz): " + ", ".join(
        f"{lo/1e3:g}k-{hi/1e3:g}k {p:.0%}" for (lo, hi), p in zip(zip(edges[:-1], edges[1:]), parts)))
    for fb in (1e6, 2e6, 4.05e6, 1e7):
        v = np.sqrt(2 * integ(lambda x: Lv(x) * x**2, 8.3e3, fb))
        print(f"  10 GHz, 8.3 kHz - {fb/1e6:5.2f} MHz: {v/1e3:6.3f} kHz rms")

    print("\n(4) FIDELITY with the measured spectrum (flicker region included), 10 GHz plane, band 8.3 kHz - 30 MHz")
    for fr in (0.5e6, 1e6, 2e6, 4e6, 10e6):
        inf = infidelity(Lv, fr, 8.3e3, 3e7)
        print(f"  fR = {fr/1e6:5.1f} MHz: 1-F = {inf:.2e} (F = {100*(1-inf):.5f}%) | 99.999% target "
              f"{'met' if inf <= 1e-5 else 'NOT met'}")
    share = infidelity(Lv, 1e6, 4e6, 1e7) / infidelity(Lv, 1e6, 8.3e3, 3e7)
    print(f"  fR = 1 MHz: share of 1-F from 4-10 MHz offsets = {share:.0%} (FN^2 share above)")

    print("\n(5) RECONCILING Gong's stated fidelity: add the measured -69 dBc reference spur at 100 MHz (inference: per sideband)")
    P, fs = 10 ** (-69 / 10), 100e6
    for fr, stated in ((1e6, 99.9994), (10e6, 99.9998)):
        base = infidelity(Lv, fr, 1e4, 3e7)
        sp = 2 * P * fs**2 / fr**2 * h_lo2(fs, fr)
        print(f"  fR = {fr/1e6:4.0f} MHz: trace {base:.2e} + spur {sp:.2e} -> F = {100*(1-base-sp):.5f}%  (Gong {stated}%)")

    print("\n(6) SENSITIVITY of the budget number (10 GHz plane, 8.3 kHz - 10 MHz)")
    fnv = lambda off: np.sqrt(2 * integ(lambda x: 10 ** ((L_db(x) + 20 * np.log10(N_DIV) + off) / 10) * x**2, 8.3e3, 1e7))
    ratio = fm / ref["FM_Hz"]
    print(f"  nominal {fnv(0)/1e3:.2f} kHz; +-0.3 dB -> {fnv(-.3)/1e3:.2f}..{fnv(.3)/1e3:.2f} kHz; "
          f"scaled by the instrument FM ratio (1/{ratio:.2f}) -> {fnv(0)/ratio/1e3:.2f} kHz")
    if os.path.exists(a.data300):
        D3 = json.load(open(a.data300))
        L3, _ = spectrum(D3)
        L3t = lambda x: 10 ** (L3(x) / 10)
        I3 = integ(L3t, 1e4, 3e7)
        fm3 = np.sqrt(2 * integ(lambda x: L3t(x) * x**2, 1e4, 3e7))
        r3 = D3["instrument_reference"]["integrated_10k_30M"]
        print(f"  300 K VALIDATION (10 kHz-30 MHz, 2.5 GHz): int noise {10*np.log10(I3):.2f} vs {r3['int_noise_dBc']} dBc, "
              f"PM {np.sqrt(2*I3)*1e6:.0f} vs {r3['PM_urad']:.0f} urad, FM {fm3:.0f} vs {r3['FM_Hz']:.0f} Hz, "
              f"jitter {np.sqrt(2*I3)/(2*np.pi*2.49991702e9)*1e15:.2f} vs {r3['jitter_fs']} fs")
        v3 = np.sqrt(2 * integ(lambda x: 10 ** ((L3(x) + 20 * np.log10(N_DIV)) / 10) * x**2, 8.3e3, 1e7))
        print(f"  300 K (Fig. 13(a)) for comparison: {v3/1e3:.2f} kHz")


if __name__ == "__main__":
    main()
