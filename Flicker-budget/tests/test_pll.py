"""Tests of the linear PLL model (src/pll.py). Run with `pytest` or `python tests/test_pll.py`.
The brick-wall test checks the implementation against the closed form; it is a unit test, not a result."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import budget as b  # noqa: E402
import pll  # noqa: E402


def test_brick_wall_limit_reproduces_closed_form():
    for al, fc in ((0.137, 4.12e6), (0.184, 0.0), (0.1, 1e6)):
        for fp in (3e5, 1e6, 3e6):
            l_in = pll.l_osc(fp, al, fc)                         # "matched" floor = oscillator at f_PLL
            got = pll.fn_rms(pll.brick_wall_h_lp(fp), al, fc, l_in, n=400001)
            want = b.fn(al, fc, fpll=fp, continuous=True)
            assert abs(got / want - 1) < 5e-3, (al, fc, fp, got, want)


def test_minus_3db_definition_and_dc_gain():
    for arch in pll.ARCH:
        for z in (0.7, 1.0, 1.6):
            h = pll.make_h_lp(arch, 1e6, z)
            assert abs(abs(h(1.0)) - 1) < 1e-6                    # unity DC gain
            assert abs(abs(h(1e6)) ** 2 - 0.5) < 1e-6             # -3 dB at the stated bandwidth


def test_high_frequency_rolloff_slopes():
    # |H|^2 ~ f^-2 for first-order and type-II (zero), ~ f^-4 for the 2nd-order low-pass
    for arch, slope in (("first-order", -2), ("type-II (PI zero)", -2), ("2nd-order low-pass", -4)):
        h = pll.make_h_lp(arch, 1e5, 1.0)
        s = np.log(abs(h(2e8)) ** 2 / abs(h(1e8)) ** 2) / np.log(2)
        assert abs(s - slope) < 0.05, (arch, s)


def test_oscillator_is_suppressed_in_band_and_untouched_far_out():
    h = pll.make_h_lp("type-II (PI zero)", 1e6)
    assert abs(1 - h(1e3)) ** 2 < 1e-5                            # deep in band: suppressed
    assert abs(abs(1 - h(1e9)) ** 2 - 1) < 1e-3                   # far out: passes through


def test_direct_spectrum_path_agrees_with_model_for_model_data():
    f = np.logspace(np.log10(1.5e5), np.log10(2.9e7), 40)
    y = 10 * np.log10(pll.l_osc(f, 0.137, 4.12e6))
    osc = pll.l_osc_from_data(f, y, 0.137, 4.12e6)
    h = pll.make_h_lp("type-II (PI zero)", 1e6)
    a = pll.fn_rms(h, 0.137, 4.12e6, 1e-13)
    c = pll.fn_rms(h, 0.137, 4.12e6, 1e-13, osc=osc)
    assert abs(a / c - 1) < 5e-3                              # interpolation error only


def test_type2_bandwidth_matches_gong_definition():
    # Gong's thesis (Ch. 5 footnote): BW = sqrt(2z^2 + 1 + sqrt((2z^2 + 1)^2 + 1)) * fn for the type-II loop
    for z in (0.7, 1.0, 1.6):
        a = 2 * z**2 + 1
        assert abs(pll.ratio_3db("type-II (PI zero)", z) / np.sqrt(a + np.sqrt(a**2 + 1)) - 1) < 1e-9
    assert pll.GONG_MATCHED == ("type-II (PI zero)", 1.6)


def test_linearity_in_input_floor():
    h = pll.make_h_lp("2nd-order low-pass", 1e6)
    o1, i1 = pll.fn_components(h, 0.137, 4.12e6, 1e-13)
    o2, i2 = pll.fn_components(h, 0.137, 4.12e6, 2e-13)
    assert abs(o1 / o2 - 1) < 1e-12 and abs(i2 / i1 - 2) < 1e-12


if __name__ == "__main__":
    names = [n for n in sorted(globals()) if n.startswith("test_")]
    for n in names:
        globals()[n]()
        print("PASS", n)
    print(f"{len(names)} tests passed")
