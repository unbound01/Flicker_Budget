#!/usr/bin/env python3
"""Digitize Fig. 7.21 (a)-(c) of J. Gong, "High-Performance Phase-Locked Loops for Quantum Computing Applications",
PhD thesis, TU Delft, 2023 (printed p. 143, PDF p. 154): measured 4.2 K phase noise of the octave-tuning DAPLL with a
dual-core triple-mode VCO (a different chip from the JSSC 2023 DAPLL), at 8.0, 10.6 and 15.4 GHz, measured after /8.

Each panel is an analyzer screenshot that prints spot-noise markers and integrated values over 10 kHz - 30 MHz;
these are stored as calibration targets. Same tracing method as digitize_gong_fig13b.py; per-panel calibration.

Usage:
    python scripts/digitize_gong_thesis_fig721.py --pdf PATH/TO/Thesis_Jiang_Gong.pdf [--out data/gong_thesis_fig721.json] [--overlay out.png]
Output is derived from a copyrighted figure: keep it out of public repos.
Calibration is for the 203-page TU Delft repository PDF, page 154 rendered at 300 dpi (2481 x 3508 px).
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import tempfile

import cv2
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import digitize_gong_fig13b as d13  # noqa: E402

X_100K, X_1M = 1108.0, 1416.0          # solid gridlines; 308 px/decade; frame edges = 10 kHz and 30 MHz
PANELS = {
    "8.0GHz": {"f_test": 1.000008437e9, "f_pll": 8.0e9, "rows": {-130: 481.0, -165: 838.0}, "plot": (431, 887),
               "ref": {"spot": {"1e4": -127.58, "1e5": -136.49, "1e6": -140.04, "1e7": -159.61},
                       "int": {"int_noise_dBc": -73.54, "PM_urad": 297.61, "FM_Hz": 1085.0, "jitter_fs": 47.365}}},
    "10.6GHz": {"f_test": 1.325013144e9, "f_pll": 10.6e9, "rows": {-140: 1177.0, -170: 1452.0}, "plot": (1041, 1497),
                "ref": {"spot": {"1e4": -126.58, "1e5": -135.21, "1e6": -137.05, "1e7": -159.43},
                        "int": {"int_noise_dBc": -72.43, "PM_urad": 338.05, "FM_Hz": 1002.0, "jitter_fs": 40.605}}},
    "15.4GHz": {"f_test": 1.925015081e9, "f_pll": 15.4e9, "rows": {-125: 1698.0, -165: 2065.0}, "plot": (1653, 2109),
                "ref": {"spot": {"1e4": -121.41, "1e5": -130.85, "1e6": -131.82, "1e7": -153.83},
                        "int": {"int_noise_dBc": -67.92, "PM_urad": 568.45, "FM_Hz": 1668.0, "jitter_fs": 46.998}}},
}


DEC = X_1M - X_100K
X_10K = X_100K - DEC
F_GRID = np.logspace(4, np.log10(3e7), 160)


def track_robust(mask, c0, c1, max_jump=10.0, reacquire=12):
    """Left-to-right; nearest cluster to a local linear prediction; after `reacquire` rejected columns the nearest
    cluster is accepted anyway (the thin thesis trace has short gaps and steep sections)."""
    path, hist, miss = {}, [], 0
    for x in range(c0, c1):
        cl = d13.clusters(mask[:, x])
        if not cl:
            continue
        if not hist:
            y = min(cl)                                  # trace starts at the top-left in all panels
        else:
            if len(hist) >= 8:
                k = np.polyfit([a for a, _ in hist[-15:]], [b for _, b in hist[-15:]], 1)
                pred, slope = np.polyval(k, x), abs(k[0])
            else:
                pred, slope = hist[-1][1], 0.0
            y = min(cl, key=lambda c: abs(c - pred))
            if abs(y - pred) > max_jump + 3.0 * slope and miss < reacquire:
                miss += 1
                continue
        miss = 0
        hist.append((x, y))
        path[x] = y
    return path


def digitize_panel(img, p):
    k, b = np.polyfit(list(p["rows"]), list(p["rows"].values()), 1)       # row = k*dB + b
    r0, r1 = p["plot"]
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV).astype(int)
    h, s, v = hsv[..., 0], hsv[..., 1], hsv[..., 2]
    m = (h >= 108) & (h <= 128) & (s > 50) & (v > 90)
    keep = np.zeros(m.shape, bool)
    keep[r0:r1, 801:1871] = True
    path = track_robust(m & keep, 802, 1870)
    xs = np.array(sorted(path), float)
    ys = np.array([path[int(x)] for x in xs])
    # 5-column running median to remove single-column outliers (marker symbols on the trace)
    ys_m = np.array([np.median(ys[max(0, i - 2):i + 3]) for i in range(len(ys))])
    f = 1e4 * 10 ** ((xs - X_10K) / DEC)
    d = (ys_m - b) / k
    sel = (F_GRID >= f.min()) & (F_GRID <= f.max())
    pts = [[float(g), float(np.interp(np.log(g), np.log(f), d))] for g in F_GRID[sel]]
    return pts, (xs, ys_m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", required=True)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "data", "gong_thesis_fig721.json"))
    ap.add_argument("--overlay", default=None)
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as td:
        subprocess.run(["pdftoppm", "-r", "300", "-f", "154", "-l", "154", "-png", a.pdf, os.path.join(td, "p")], check=True)
        img = cv2.imread(glob.glob(os.path.join(td, "p-*.png"))[0])
    out = {"source": "Gong PhD thesis (TU Delft 2023) Fig. 7.21, 4.2 K, octave DAPLL, measured after /8",
           "units": "offset [Hz], SSB PN [dBc/Hz] at the test output f_test = f_pll/8; add 20log10(8) = 18.06 dB for f_pll",
           "panels": {}}
    ov = img.copy() if a.overlay else None
    for name, p in PANELS.items():
        pts, (xs, ys) = digitize_panel(img, p)
        t = np.array(pts)
        print(f"{name}: {len(pts)} points, {t[0,0]/1e3:.1f} kHz - {t[-1,0]/1e6:.1f} MHz")
        for k, v in p["ref"]["spot"].items():
            print(f"   marker {float(k):7.0e}: digitized {np.interp(np.log(float(k)), np.log(t[:,0]), t[:,1]):8.2f}  instrument {v}")
        out["panels"][name] = {"f_test": p["f_test"], "f_pll": p["f_pll"], "instrument_reference": p["ref"], "trace": pts}
        if ov is not None:
            for x, y in zip(xs[::3], ys[::3]):
                cv2.circle(ov, (int(x), int(y)), 2, (0, 0, 255), -1)
    json.dump(out, open(a.out, "w"), indent=1)
    print("wrote", a.out)
    if ov is not None:
        cv2.imwrite(a.overlay, ov[350:2210, 780:1890])


if __name__ == "__main__":
    main()
