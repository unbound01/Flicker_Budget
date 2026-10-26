#!/usr/bin/env python3
"""Three criteria on the same LO spectra: is 1.9 kHz the right number for any operation?

Spectra (single-sideband L(f), 1/Hz):
  G4   Gong JSSC 2023 Fig. 13(b), measured PLL output at 4.2 K, referred to the ~10 GHz LO (+12.04 dB)
  G300 Gong Fig. 13(a), same at 300 K
  P-FR Patra JSSC 2018 4 K oscillator, two-term fit alpha/f^2 (1 + fc/f), free running (6.3 GHz)
  P-L1, P-L4  the same oscillator in the Gong-matched loop (type-II, zeta = 1.6), BW 1 / 4 MHz, white-FM floor alpha/f3^2
  T-8.0, T-10.6, T-15.4  second chip: Gong thesis Fig. 7.21, octave DAPLL, measured at 4.2 K, referred to f_pll (+18.06 dB)
  X-8.6  third chip: Xiao et al. RFIC 2025, 28 nm fractional-N CP-PLL, integer mode, 8.6 GHz, 4 K
  U-5.7 ... U-6.9  fourth chip: Peng et al. RFIT 2024, 65 nm sub-sampling PLL, 4 K, at the plotted output frequency

Criteria (S_f = 2 f^2 L is the one-sided frequency-noise PSD, Hz^2/Hz):
  brick   Patra:  sigma^2 = int_{8.3k}^{10M} S_f df                           compared with 1.9 kHz
  gate    Gong Eqs. 1-2 (pi rotation):  1-F = int S_f |H_LO|^2 / fR^2 df       compared with 1e-5
  idle    free evolution (Ramsey) of length tau: accumulated phase variance
             var_phi = (2 pi tau)^2 int S_f sinc^2(pi f tau) df,  sinc(x) = sin(x)/x
          reported as the effective quasi-static FN seen by the qubit, sigma_eff(tau) = sqrt(var_phi)/(2 pi tau),
          which is directly comparable with the qubit's own quasi-static 1.9 kHz (T2* = 1/(sqrt2 pi sigma) = 118 us).

Below the first measured point (10.5 kHz) the Gong spectra need an assumption; three scenarios are run:
  A  continue the measured local slope (about -9 dB/dec, i.e. close to flicker PM, which Gong attributes to the INBUF)
  B  white FM below 10 kHz (-20 dB/dec): S_f flat
  C  pessimistic: A down to 1 kHz, then flicker FM (-30 dB/dec, S_f ~ 1/f) down to F_MIN
Lowest frequency F_MIN = 1 Hz (a ~1 s averaging window); slower drift is treated as calibrated out.

Usage: python scripts/filter_criteria.py [--g4 data/gong_fig13b.json] [--g300 data/gong_fig13a.json] [--patra data/patra_digitized.json]
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
import pll  # noqa: E402
from a0_gong_measured import h_lo2  # noqa: E402

F_MIN, F_MAX = 1.0, 1e8                 # integration band [Hz]; change with set_band() or --f-min/--f-max
GRID = np.logspace(np.log10(F_MIN), np.log10(F_MAX), 600001)
LNF = np.log(GRID)


def set_band(f_min=1.0, f_max=1e8, n=600001):
    """Set the integration band used by every criterion (gate, idle, and the spectra built on GRID)."""
    global F_MIN, F_MAX, GRID, LNF
    F_MIN, F_MAX = float(f_min), float(f_max)
    GRID = np.logspace(np.log10(F_MIN), np.log10(F_MAX), n)
    LNF = np.log(GRID)
DB_DIV4 = 20 * np.log10(4)


def integ(y):
    """int y(f) df on GRID (log grid: df = f dlnf)."""
    return np.trapezoid(y * GRID, LNF)


def sinc2_avg(f, tau):
    """sinc^2(pi f tau), replaced by its period average 1/(2 (pi f tau)^2) where the grid cannot resolve it."""
    x = np.pi * f * tau
    out = np.where(x < 1e-8, 1.0, (np.sin(x) / np.where(x == 0, 1, x)) ** 2)
    return np.where(x > 200, 0.5 / x**2, out)


# ---- criteria -----------------------------------------------------------------------------------------
def brick(L):
    m = (GRID >= 8.3e3) & (GRID <= 1e7)
    return np.sqrt(np.trapezoid((2 * GRID**2 * L * GRID)[m], LNF[m]))


def gate(L, fr):
    return integ(2 * GRID**2 * L * h_lo2(GRID, fr) / fr**2)


def idle_sigma_eff(L, tau):
    return np.sqrt(integ(2 * GRID**2 * L * sinc2_avg(GRID, tau)))


def idle_var(L, tau):
    """Ramsey (free-evolution) phase variance [rad^2], Eq. (9) of the paper."""
    return (2 * np.pi * tau) ** 2 * idle_sigma_eff(L, tau) ** 2


def echo_var(L, tau):
    """Hahn-echo phase variance [rad^2] for total duration tau (pi pulse at tau/2, ideal pulses).
    Phase-error filter |1 - 2 e^{j pi f tau} + e^{j 2 pi f tau}|^2 = 16 sin^4(pi f tau / 2); S_phi = 2 L.
    Above ~200 rad the weight is replaced by its cycle average, 6."""
    x = np.pi * GRID * tau / 2
    return integ(2 * L * np.where(x > 200, 6.0, 16 * np.sin(x) ** 4))


def idle_infidelity(var, worst_case=False):
    """Gaussian phase error: average infidelity (1-e^{-v/2})/3, worst case (equatorial state) (1-e^{-v/2})/2."""
    return (1 - np.exp(-var / 2)) / (2 if worst_case else 3)


# ---- spectra ------------------------------------------------------------------------------------------
def gong_spectrum(path, scenario, hi="slope"):
    return trace_spectrum(np.array(json.load(open(path))["trace"]), DB_DIV4, scenario, hi)


def trace_spectrum(t, db_offset, scenario="A", hi="slope"):
    """Digitized trace [[f, dBc/Hz], ...] -> L(GRID) (linear), referred by db_offset, with the low-offset scenario.
    hi = "slope": beyond the last point continue the end slope (default); "flat": hold the last measured level."""
    f, d = t[:, 0], t[:, 1] + db_offset
    lf = np.log10(f)
    k_lo = np.polyfit(lf[:6], d[:6], 1)
    k_hi = np.polyfit(lf[-6:], d[-6:], 1)
    lg = np.log10(GRID)
    out = np.interp(lg, lf, d)
    above = GRID > f[-1]
    out[above] = np.polyval(k_hi, lg[above]) if hi == "slope" else d[-1]
    lo = GRID < f[0]
    d0 = np.polyval(k_lo, np.log10(f[0]))
    if scenario == "A":
        out[lo] = np.polyval(k_lo, lg[lo])
    elif scenario == "B":
        out[lo] = d0 - 20 * (lg[lo] - np.log10(f[0]))
    elif scenario == "C":
        out[lo] = np.polyval(k_lo, lg[lo])
        m = GRID < 1e3
        d1k = np.polyval(k_lo, 3.0)
        out[m] = d1k - 30 * (lg[m] - 3.0)
    return 10 ** (out / 10), k_lo[0]


def patra_free(alpha, fc):
    return pll.l_osc(GRID, alpha, fc)


def patra_in_loop(alpha, fc, bw):
    arch, zeta = pll.GONG_MATCHED
    h = pll.make_h_lp(arch, bw, zeta)(GRID)
    return np.abs(h) ** 2 * (alpha / bw**2) + np.abs(1 - h) ** 2 * pll.l_osc(GRID, alpha, fc)


# ---- self-checks (unit tests of the filter code, not results) -------------------------------------------
def self_checks():
    # 1) quasi-static limit: a narrow band of FN at very low f is seen fully by the idle filter
    L = np.where((GRID > 1.0) & (GRID < 2.0), 1.0, 0.0) / (2 * GRID**2)   # S_f = 1 Hz^2/Hz on 1-2 Hz -> sigma^2 = 1 Hz^2
    assert abs(idle_sigma_eff(L, 1e-6) - 1.0) < 1e-3
    # 2) qubit's own number: sigma = 1.9 kHz quasi-static gives coherence exp(-var/2) = exp(-(tau/T2*)^2) at T2*
    t2 = 1 / (np.sqrt(2) * np.pi * 1.9e3)
    assert abs(t2 * 1e6 - 118.5) < 0.5
    # 3) white PM far above 1/tau: var_phi -> 2 * (2 int L df) = 2 * sigma_phi^2 (independent start/end phase errors)
    L = np.where((GRID > 1e6) & (GRID < 1e7), 1e-14, 0.0)
    var = (2 * np.pi * 1e-4) ** 2 * idle_sigma_eff(L, 1e-4) ** 2
    assert abs(var / (2 * 2 * 1e-14 * 9e6) - 1) < 0.02
    # 4) gate filter at DC is 1 and the white-FM pi-pulse result matches the repo's report value
    assert abs(h_lo2(np.array([1e-3]), 1e6)[0] - 1) < 1e-9
    L = 0.184 / GRID**2
    assert abs(gate(L, 1e6) / 9.08e-7 - 1) < 0.01
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--g4", default=os.path.join(HERE, "..", "data", "gong_fig13b.json"))
    ap.add_argument("--g300", default=os.path.join(HERE, "..", "data", "gong_fig13a.json"))
    ap.add_argument("--patra", default=os.path.join(HERE, "..", "data", "patra_digitized.json"))
    ap.add_argument("--thesis", default=os.path.join(HERE, "..", "data", "gong_thesis_fig721.json"))
    ap.add_argument("--xiao", default=os.path.join(HERE, "..", "data", "xiao_rfic2025_fig7.json"))
    ap.add_argument("--peng", default=os.path.join(HERE, "..", "data", "peng_rfit2024_fig6.json"))
    ap.add_argument("--f-min", type=float, default=1.0, help="lower integration limit [Hz]")
    ap.add_argument("--f-max", type=float, default=1e8, help="upper integration limit [Hz]")
    a = ap.parse_args()
    set_band(a.f_min, a.f_max)
    print(f"integration band for gate and idle: {F_MIN:g} Hz - {F_MAX:g} Hz (brick: 8.3 kHz - 10 MHz)")
    print("self-checks:", "PASS" if self_checks() else "FAIL", "(unit tests of the filter code, not results)")

    alpha, fc = 0.1372, 4.122e6                 # primary 4 K fit (docs/identifiability.md); refit if data present
    if os.path.exists(a.patra):
        sys.path.insert(0, HERE)
        from pll_comparison import fit_primary
        A = np.array(json.load(open(a.patra))["fig22_PN_dBc_per_Hz_vs_offset_Hz"]["4K"])
        alpha, fc = fit_primary(A)
    spectra = {}
    for sc in "ABC":
        L, k = gong_spectrum(a.g4, sc)
        spectra[f"G4-{sc}"] = L
    spectra["G300-A"], _ = gong_spectrum(a.g300, "A")
    if os.path.exists(a.thesis):
        for name, pnl in json.load(open(a.thesis))["panels"].items():
            spectra["T-" + name.replace("GHz", "")], _ = trace_spectrum(np.array(pnl["trace"]), 20 * np.log10(8), "A")
    if os.path.exists(a.xiao):                         # independent chip 3: Xiao et al. RFIC 2025, 8.6 GHz, no divider
        spectra["X-8.6"], _ = trace_spectrum(np.array(json.load(open(a.xiao))["trace"]), 0.0, "A")
    if os.path.exists(a.peng):                         # independent chip 4: Peng et al. RFIT 2024, plotted at the output
        for name, pnl in json.load(open(a.peng))["panels"].items():
            spectra["U-" + name.replace("GHz", "")], _ = trace_spectrum(np.array(pnl["trace"]), 0.0, "A")
    spectra["P-FR"] = patra_free(alpha, fc)
    spectra["P-L1"] = patra_in_loop(alpha, fc, 1e6)
    spectra["P-L4"] = patra_in_loop(alpha, fc, 4e6)
    print(f"Gong 4.2 K low-offset slope used in scenario A: {k:.1f} dB/dec; Patra fit alpha={alpha:.4f}, fc={fc/1e6:.2f} MHz")

    taus = (1e-6, 1e-5, 1.2e-4)
    print("\n spectrum | brick FN 8.3k-10M | gate 1-F fR=1 / 10 MHz | idle sigma_eff [Hz] at tau = 1 / 10 / 120 us | idle at 120 us / 1.9 kHz")
    rows = {}
    for name, L in spectra.items():
        b = brick(L)
        g1, g10 = gate(L, 1e6), gate(L, 1e7)
        idl = [idle_sigma_eff(L, t) for t in taus]
        rows[name] = (b, g1, g10, idl)
        print(f" {name:7s} | {b/1e3:7.2f} kHz ({b/1.9e3:4.1f}x) | {g1:.1e} / {g10:.1e} | "
              + " / ".join(f"{x:9.1f}" for x in idl) + f" | {idl[2]/1.9e3:.1e}")

    print("\nIdle as accumulated phase variance var_phi(tau) [rad^2], against the qubit's own quasi-static noise"
          " (2 pi 1.9 kHz tau)^2; and 2*sigma_phi^2 (integrated PM, start + end) as the locked-LO asymptote:")
    q = lambda t: (2 * np.pi * 1.9e3 * t) ** 2
    print("  tau [us]        " + "  ".join(f"{t*1e6:>16.1f}" for t in (1e-7, 1e-6, 1e-5, 1.2e-4)) + "   2*sigma_phi^2")
    print("  qubit's own     " + "  ".join(f"{q(t):16.2e}" for t in (1e-7, 1e-6, 1e-5, 1.2e-4)))
    for name, L in spectra.items():
        v = [(2 * np.pi * t) ** 2 * idle_sigma_eff(L, t) ** 2 for t in (1e-7, 1e-6, 1e-5, 1.2e-4)]
        print(f"  {name:15s} " + "  ".join(f"{x:9.2e} ({x/q(t):5.0e})" for x, t in zip(v, (1e-7, 1e-6, 1e-5, 1.2e-4)))
              + f"   {2 * 2 * integ(L):.2e}")

    print("\nWhere each criterion looks (Gong 4.2 K, scenario A): share of the criterion from offset bands")
    L = spectra["G4-A"]
    bands = [(F_MIN, 1e3), (1e3, 1e4), (1e4, 1e5), (1e5, 1e6), (1e6, 4e6), (4e6, 1e7), (1e7, F_MAX)]

    def share(w):
        tot = integ(w)
        return [np.trapezoid((w * GRID)[(GRID >= lo) & (GRID <= hi)], LNF[(GRID >= lo) & (GRID <= hi)]) / tot for lo, hi in bands]
    sf = 2 * GRID**2 * L
    crit = {"brick (8.3k-10M)": np.where((GRID >= 8.3e3) & (GRID <= 1e7), sf, 0.0),
            "gate fR = 1 MHz": sf * h_lo2(GRID, 1e6),
            "idle tau = 10 us": sf * sinc2_avg(GRID, 1e-5),
            "idle tau = 120 us": sf * sinc2_avg(GRID, 1.2e-4)}
    print("  bands: " + " | ".join(f"{lo:.0e}-{hi:.0e}" for lo, hi in bands))
    for k, w in crit.items():
        print(f"  {k:18s}: " + " | ".join(f"{s:7.1%}" for s in share(w)))

    print("\nIdle sigma_eff(tau) curve, Gong 4.2 K (A / B / C) and Patra free-running, tau in us:")
    for t in (0.1e-6, 1e-6, 1e-5, 1e-4, 1.2e-4, 1e-3):
        print(f"  tau {t*1e6:7.1f}: " + "  ".join(f"{n} {idle_sigma_eff(spectra[n], t):9.1f} Hz" for n in ("G4-A", "G4-B", "G4-C", "P-FR", "P-L4")))


if __name__ == "__main__":
    main()
