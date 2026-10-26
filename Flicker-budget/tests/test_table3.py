"""Regression test: scripts/table3.py reproduces the published Table III anchors (skipped without local data)."""
import json
import os
import subprocess
import sys
import tempfile

import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
NEEDED = ["gong_fig13b", "gong_fig13a", "gong_thesis_fig721", "xiao_rfic2025_fig7", "peng_rfit2024_fig6", "patra_digitized"]
HAVE = all(os.path.exists(os.path.join(ROOT, "data", n + ".json")) for n in NEEDED)


@pytest.mark.skipif(not HAVE, reason="digitized data not present (regenerate with the digitize_* scripts)")
def test_table3_anchors():
    with tempfile.TemporaryDirectory() as d:
        subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "table3.py"), "--out", d], check=True, capture_output=True)
        N = json.load(open(os.path.join(d, "numbers.json")))
        rows = open(os.path.join(d, "table3_rows.tex")).read()
    p = N["points"]
    assert abs(p["G7-10.0"]["own"]["brick"] - 8.86e3) < 50
    assert abs(p["G7-10.0"]["own"]["gate"]["1"] / 5.5e-6 - 1) < 0.03
    assert abs(p["X12-8.6"]["own"]["gate"]["1"] / 7.9e-7 - 1) < 0.03
    assert sorted(N["fail_echo_avg_own"]) == ["G10-15.4", "P13-5.7", "P13-6.0", "P13-6.3", "P13-6.9"]
    assert N["discordant_pairs"] == [8, 36]
    assert abs(N["free_running"]["patra_6.3"]["ramsey_avg"] - 0.309) < 0.005
    assert abs(N["jitter_spec"]["echo_fs_10GHz"] - 50.3) < 0.1
    assert rows.count("\\\\") == 9
    lp = N["patra_in_loop"]
    assert abs(lp["1e+06"]["white_FM"]["gate1"] / 2.0e-6 - 1) < 0.05 and abs(lp["4e+06"]["flat_-122"]["gate1"] / 6.8e-6 - 1) < 0.03
    w = N["weights_gong7"]
    assert abs(w["brick_4_10MHz"] - 0.77) < 0.01 and abs(w["gate1_1_4MHz"] - 0.66) < 0.01
    assert abs(w["fidelity_fR1"] - 0.9999945) < 2e-7
    fr = N["free_running"]
    assert fr["gong2022_6.3"]["ramsey_tracked_avg"] > 40 * 1e-5 and fr["gong2022_6.3"]["echo_avg"] > 40 * 1e-5   # "at least 40-fold"
    assert abs(lp["1e+06"]["white_FM"]["ramsey_avg"] / 1.5e-6 - 1) < 0.05


if __name__ == "__main__":
    test_table3_anchors() if HAVE else print("skipped: no data")
    print("ok")
