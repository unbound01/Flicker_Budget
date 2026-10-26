"""Regression tests. Run with `pytest` or simply `python tests/test_budget.py`.

Note: the closed-form-vs-quad test is a unit test of the algebra, not a scientific result.
"""
import json
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import budget as b  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

DATA_PATH = os.path.join(HERE, "..", "data", "patra_digitized.json")
DATA = json.load(open(DATA_PATH)) if os.path.exists(DATA_PATH) else None   # generated locally; not distributed


def need_data():
    if DATA is None:
        raise unittest.SkipTest("data/patra_digitized.json not generated (see data/README.md)")


def approx(x, y, tol):
    return abs(x - y) <= tol


# ---- baseline (Patra Eqs. 7-9) ---------------------------------------------------------
def test_alpha_max_and_phase_noise():
    a = b.alpha_max()
    assert approx(a, 0.1842, 5e-4)
    assert approx(b.dbc(a / b.FPLL**2), -116.9, 0.1)     # in-band (Patra quotes -115: open item)
    assert approx(b.dbc(a / b.FB**2), -147.3, 0.1)       # at 10 MHz (Patra: -147)
    assert approx(b.fn(0.2) / 1e3, 1.98, 0.01)           # Patra's alpha = 0.2 slightly exceeds 1.9 kHz


def test_bandwidth_must_lie_inside_integration_band():
    for bad in (5e3, 2e7):
        try:
            b.fn(0.1, 1e5, fpll=bad)
        except ValueError:
            continue
        raise AssertionError(f"no ValueError for fpll={bad:g}")


def test_closed_form_matches_numerical_integration():
    assert b.unit_test() < 1e-3


def test_corner_tolerance_values():
    assert approx(b.fc_max(0.17) / 1e3, 233.2, 1)
    assert approx(b.fc_max(0.15) / 1e3, 636.9, 1)
    assert approx(b.fc_max(0.10) / 1e3, 2353, 5)
    assert approx(b.fc_max(0.17, continuous=True) / 1e3, 212.9, 1)


def test_flicker_budget_depends_only_on_alpha_and_K():
    # same K = alpha*fc, different alpha: the flicker term must be identical
    K = 5e4
    base = lambda al: b.fn2(al, K / al) - b.fn2(al, 0.0)
    assert approx(base(0.05), base(0.15), 1e-6 * base(0.05))


def test_fn_example_and_breakdown():
    a = b.alpha_max()
    assert approx(b.fn(a, 500e3) / 1e3, 2.063, 2e-3)


def test_fb_table():
    got = [b.alpha_max(fb=f) for f in (1e6, 2e6, 5e6, 10e6)]
    assert np.allclose(got, [2.26, 1.00, 0.38, 0.18], atol=0.006)


# ---- unit / convention checks ----------------------------------------------------------
def test_van_dijk_table_I():
    L = 10 ** (-10.6)                                    # L(1 MHz) = -106 dBc/Hz
    sigma = np.sqrt(2 * L * 1e12 * 2.5e6)                # FN^2 = 2 L f^2 ENBW
    assert approx(sigma / 1e3, 11.2, 0.1)                # table: 11 kHz
    assert approx((1.9e3 / 1e6) ** 2, 3.6e-6, 0.05e-6)   # 1.9 kHz at fR = 1 MHz


def test_ball_filter_matches_van_dijk_for_white_fm():
    from scipy.integrate import quad

    def prim(z):                                         # Ball SI Eq. 112
        d = (np.pi**2 - z**2) ** 2
        return (2 * np.pi**2 * z**2 + 2 * z**4) * (1 + np.cos(z)) / d

    alpha, fR = 0.184, 1e6
    tau = 1 / (2 * fR)
    chi = alpha * quad(lambda f: float(prim(2 * np.pi * f * tau)) / f**2, 1e2, 1e12,
                       points=[1e5, 1e6, 1e7, 1e8], limit=2000)[0]
    chi += alpha * (4 / np.pi**2) * (2 * np.pi * tau) ** 2 * 1e2
    assert approx(0.5 * chi / (alpha * np.pi**2 / (2 * fR)), 1.0, 0.01)


# ---- fit to digitized Patra data ---------------------------------------------------------
def _fit_4k(lo=2e5, hi=1e7, shift=0.0):
    a = np.array(DATA["fig22_PN_dBc_per_Hz_vs_offset_Hz"]["4K"])
    a = a[(a[:, 0] >= lo) & (a[:, 0] <= hi * 1.001)]
    f, y = a[:, 0], a[:, 1] + shift
    model = lambda p: 10 * np.log10(np.exp(p[0]) / f**2 * (1 + np.exp(p[1]) / f)) - y
    r = least_squares(model, np.log((0.3, 1e6)))
    return (*np.exp(r.x), np.sqrt(np.mean(r.fun**2)))


def test_4k_fit_reproduces_3p4_kHz_and_fails_1MHz_pll_claim():
    need_data()
    al, fc, rms = _fit_4k()
    assert rms < 1.0
    assert 3.5e6 < fc < 7e6                               # Fig. 24: 5.7 MHz at 4 K
    free = np.sqrt(b.fn2_free_running(al, fc)) / 1e3
    assert 3.0 < free < 3.6                               # Patra: 3.4 kHz
    for cont in (False, True):
        assert b.fn(al, fc, fpll=1e6, continuous=cont) / 1e3 > 1.9   # claim not supported


def test_fit_sensitivity_band():
    need_data()
    outs = []
    for lo, hi, sh in ((2e5, 1e7, 0), (1.5e5, 2.9e7, 0), (3e5, 1e7, 0), (2e5, 1e7, .5), (2e5, 1e7, -.5)):
        al, fc, _ = _fit_4k(lo, hi, sh)
        outs.append(b.fn(al, fc, fpll=1e6) / 1e3)
    assert 2.0 < min(outs) and max(outs) < 2.6


def test_range_grid_pll_conclusion_is_robust():
    need_data()
    A = np.array(DATA["fig22_PN_dBc_per_Hz_vs_offset_Hz"]["4K"])
    fn1, K = [], []
    for lo in (1.5e5, 2e5, 3e5, 5e5):
        for hi in (5e6, 1e7, 2e7, 2.9e7):
            m = (A[:, 0] >= lo) & (A[:, 0] <= hi * 1.001)
            f, y = A[m, 0], A[m, 1]
            r = least_squares(lambda p: 10 * np.log10(np.exp(p[0]) / f**2 * (1 + np.exp(p[1]) / f)) - y,
                              np.log((0.3, 1e6)))
            al, fc = np.exp(r.x)
            fn1.append(b.fn(al, fc, fpll=1e6) / 1e3)
            K.append(al * fc)
    assert min(fn1) > 1.9                      # never reaches the budget with a 1 MHz PLL
    assert max(K) / min(K) < 2.0               # K is identified far better than alpha or fc individually


def test_digitized_data_shape():
    need_data()
    assert len(DATA["fig24_corner_Hz_vs_T_K"]) == 11
    assert all(len(v) == 8 for v in DATA["fig25_PN_dBc_per_Hz_vs_GHz"].values())


if __name__ == "__main__":                                # works without pytest
    names = [n for n in sorted(globals()) if n.startswith("test_")]
    ok = skipped = 0
    for n in names:
        try:
            globals()[n]()
            print("PASS", n)
            ok += 1
        except unittest.SkipTest as e:
            print("SKIP", n, "-", e)
            skipped += 1
    print(f"{ok} passed, {skipped} skipped")
