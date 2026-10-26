#!/usr/bin/env python3
# Usage: python scripts/alternative_integrator.py data
# Alternative integrator: recomputes brick, gate and idle with integration code independent of src/ and
# filter_criteria.py. It uses the same digitized traces, so it checks the numerics, not the digitization.
"""Alternative-integrator recomputation of the headline numbers. Uses NO repo function: own interpolation, own integration
(scipy.integrate.quad on piecewise segments + dense linear-grid Simpson for the oscillatory idle filter), own
filters written from the papers' equations. Two extrapolation choices are compared."""
import json, sys
import numpy as np
from scipy.integrate import quad, simpson

D = sys.argv[1]
J = lambda f: json.load(open(f"{D}/{f}.json"))


def make_L(trace, offset_db, extrap):
    t = np.array(trace)
    lf, d = np.log10(t[:, 0]), t[:, 1] + offset_db
    s_lo = (d[3] - d[0]) / (lf[3] - lf[0])
    s_hi = (d[-1] - d[-4]) / (lf[-1] - lf[-4])
    def L(f):
        x = np.log10(f)
        if x < lf[0]:
            y = d[0] + (s_lo * (x - lf[0]) if extrap == "slope" else 0.0)          # 'flat' = hold the first value
        elif x > lf[-1]:
            y = d[-1] + (s_hi if extrap == "slope" else -20.0) * (x - lf[-1])       # 'flat' = -20 dB/dec beyond
        else:
            y = np.interp(x, lf, d)
        return 10 ** (y / 10)
    return L, t[0, 0], t[-1, 0]


def qlog(fun, a, b):
    """int_a^b fun(f) df via quad in ln f, split per half-decade."""
    edges = np.logspace(np.log10(a), np.log10(b), int(2 * np.log10(b / a)) + 2)
    return sum(quad(lambda u: fun(np.exp(u)) * np.exp(u), np.log(x0), np.log(x1), limit=200, epsrel=1e-8)[0]
               for x0, x1 in zip(edges[:-1], edges[1:]))


def H2(f, fr):                                    # Gong JSSC 2023 Eq. 2, written from the paper
    u = f / fr
    if abs(u - 1) < 1e-6:
        u = 1 + 1e-6
    return (1 + np.cos(np.pi * u)) * (1 + u * u) / (2 * (1 - u * u) ** 2)


def brick(L):
    return np.sqrt(qlog(lambda f: 2 * f * f * L(f), 8.3e3, 1e7))


def gate(L, fr):                                  # Gong Eq. 1, integrated 1 Hz .. 1 GHz in segments
    return qlog(lambda f: 2 * L(f) * f * f / fr ** 2 * H2(f, fr), 1.0, 1e9)


def idle_var(L, tau):                             # (2 pi tau)^2 int S_f sinc^2(pi f tau) df, S_f = 2 f^2 L
    lo = np.logspace(0, np.log10(20 / tau), 20001)                          # resolved part on a log grid
    yl = np.array([2 * f * f * L(f) for f in lo]) * np.sinc(lo * tau) ** 2  # np.sinc(x) = sin(pi x)/(pi x)
    part1 = simpson(yl * lo, x=np.log(lo))
    part2 = qlog(lambda f: 2 * f * f * L(f) * 0.5 / (np.pi * f * tau) ** 2, 20 / tau, 1e9)   # cycle-averaged tail
    return (2 * np.pi * tau) ** 2 * (part1 + part2)


spec = {}
g = J("gong_fig13b"); spec["Gong JSSC 10 GHz 4.2K"] = (g["trace"], 20 * np.log10(4))
t = J("gong_thesis_fig721")
for n, p in t["panels"].items():
    spec[f"Gong thesis {n} 4.2K"] = (p["trace"], 20 * np.log10(8))
x = J("xiao_rfic2025_fig7"); spec["Xiao 8.6 GHz 4K"] = (x["trace"], 0.0)
pg = J("peng_rfit2024_fig6")
for n, p in pg["panels"].items():
    spec[f"Peng {n} 4K"] = (p["trace"], 0.0)

print(f"{'spectrum':26s} {'extrap':6s} {'brick kHz':>9s} {'gate fR=1M':>11s} {'gate fR=4M':>11s} {'idle 120us':>11s}")
for name, (tr, off) in spec.items():
    for ex in ("slope", "flat"):
        L, f0, f1 = make_L(tr, off, ex)
        print(f"{name:26s} {ex:6s} {brick(L)/1e3:9.2f} {gate(L,1e6):11.2e} {gate(L,4e6):11.2e} {idle_var(L,1.2e-4):11.2e}")

# convention check: the Gong 300 K analyzer box from first principles
a = J("gong_fig13a"); L, _, _ = make_L(a["trace"], 0.0, "slope")
I = qlog(L, 1e4, 3e7)
print("\n300 K box check: int noise %.2f dBc (-61.55), PM %.0f urad (1180), FM %.0f Hz (5949)" %
      (10 * np.log10(I), np.sqrt(2 * I) * 1e6, np.sqrt(qlog(lambda f: 2 * f * f * L(f), 1e4, 3e7))))
# T2* sanity
print("T2* from 1.9 kHz quasi-static: %.1f us" % (1e6 / (np.sqrt(2) * np.pi * 1.9e3)))
# gate filter sanity: white FM alpha=0.184 at fR = 1 MHz (report: 9.08e-7)
print("white-FM gate check: %.3e (9.08e-7)" % gate(lambda f: 0.184 / f ** 2, 1e6))
