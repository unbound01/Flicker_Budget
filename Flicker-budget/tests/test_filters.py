"""Unit tests of the filter-function code in scripts/filter_criteria.py (no data needed). Not results."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import filter_criteria as F  # noqa: E402


def test_echo_refocuses_free_running_but_not_locked():
    tau = 1.2e-4
    fr = F.patra_free(0.041, 1.3e6)                    # free-running: echo removes most of the variance
    assert F.echo_var(fr, tau) < 0.15 * F.idle_var(fr, tau)
    lk = F.patra_in_loop(0.137, 4.12e6, 1e6)           # locked: high-offset weight 6 vs 2, echo is ~3x worse
    r = F.echo_var(lk, tau) / F.idle_var(lk, tau)
    assert 2.5 < r < 3.5
    v = 1e-3
    assert abs(F.idle_infidelity(v, True) / F.idle_infidelity(v) - 1.5) < 1e-9


def test_filter_self_checks():
    # quasi-static limit, T2* from 1.9 kHz, white-PM limit (2 sigma_phi^2), gate filter DC gain and white-FM value
    assert F.self_checks()


def test_locked_lo_idle_error_is_tau_independent_phase_jitter():
    # a PLL-like spectrum (flat in-band PN, f^-2 beyond 1 MHz): idle phase variance saturates at 2 * 2 int L df
    L = np.where(F.GRID < 1e6, 1e-13, 1e-13 * (1e6 / F.GRID) ** 2)
    asym = 2 * 2 * F.integ(L)
    for tau in (1e-5, 1.2e-4):
        v = (2 * np.pi * tau) ** 2 * F.idle_sigma_eff(L, tau) ** 2
        assert abs(v / asym - 1) < 0.05


def test_free_running_flicker_fm_idle_noise_grows_with_tau():
    L = F.pll.l_osc(F.GRID, 0.137, 4.12e6)
    v1 = (2 * np.pi * 1e-5) ** 2 * F.idle_sigma_eff(L, 1e-5) ** 2
    v2 = (2 * np.pi * 1e-4) ** 2 * F.idle_sigma_eff(L, 1e-4) ** 2
    assert v2 / v1 > 50          # ~tau^2 (quasi-static FN), unlike a locked LO


if __name__ == "__main__":
    names = [n for n in sorted(globals()) if n.startswith("test_")]
    for n in names:
        globals()[n]()
        print("PASS", n)
    print(f"{len(names)} tests passed")
