#!/usr/bin/env python3
"""Digitize Fig. 13(b) of J. Gong et al., "A Cryo-CMOS PLL for Quantum Computing Applications", IEEE JSSC 58(5), 2023
(journal p. 1371): phase-noise analyzer screenshot, PLL locked at ~10 GHz, 4.2 K, measured after the on-chip divide-by-4
(signal frequency 2.500018 GHz).

The screenshot also prints the instrument's own results, used here as calibration targets:
  spot noise  10 kHz -127.03, 100 kHz -132.74, 1 MHz -133.87, 10 MHz -149.01 dBc/Hz
  integrated 10 kHz - 30 MHz: -67.80 dBc, PM 575.92 urad, FM 2.308 kHz, jitter 36.663 fs

Usage:
    python scripts/digitize_gong_fig13b.py --pdf PATH/TO/Gong_JSSC2023.pdf [--out data/gong_fig13b.json] [--overlay out.png]

You must supply your own copy of the paper. The output is derived from a copyrighted figure: keep it out of public repos.
Calibration is for the 14-page IEEE Xplore PDF, page 10 rendered at 300 dpi (2550 x 3300 px).
"""
import argparse
import glob
import json
import os
import subprocess
import tempfile

import cv2
import numpy as np

# ---- calibration (page 10 @ 300 dpi) ------------------------------------------------------------
X_100K, X_1M, X_10M = 512.0, 806.0, 1100.0        # solid vertical gridlines; 294 px/decade
DEC = (X_10M - X_100K) / 2.0
X_10K = X_100K - DEC                              # = 218, plot frame at 216-217
GRID_ROWS = {-130: 993.0, -135: 1038.0, -145: 1126.0, -155: 1212.0, -165: 1300.0}   # horizontal gridlines
_k, _b = np.polyfit(list(GRID_ROWS.keys()), list(GRID_ROWS.values()), 1)            # row = _k*dB + _b
PLOT = (952, 1386, 219, 1240)                     # rows, cols inside the trace area (below the header bar)
EXCLUDE = [(1048, 1100, 300, 440)]                # "4.2K" label box (same blue as the trace)
F_GRID = np.logspace(np.log10(1.0e4), np.log10(3.0e7), 160)


def f_of(x):
    return 1e4 * 10 ** ((np.asarray(x, float) - X_10K) / DEC)


def x_of(f):
    return X_10K + np.log10(np.asarray(f) / 1e4) * DEC


def db_of(y):
    return (np.asarray(y, float) - _b) / _k


def y_of(db):
    return _k * np.asarray(db) + _b


def clusters(col):
    idx = np.nonzero(col)[0]
    if len(idx) == 0:
        return []
    out, st, pv = [], idx[0], idx[0]
    for i in idx[1:]:
        if i - pv > 2:
            out.append((st + pv) / 2)
            st = i
        pv = i
    out.append((st + pv) / 2)
    return out


def track(mask, x_start, x_stop, max_jump=6.0):
    """Follow the trace from the right edge leftwards; rejects jumps (marker labels, the narrow spur)."""
    hist, path = [], {}
    for x in range(x_start, x_stop, -1):
        cl = clusters(mask[:, x])
        if not cl:
            continue
        if len(hist) >= 10:
            px = np.array([a for a, _ in hist[-20:]], float)
            py = np.array([b for _, b in hist[-20:]], float)
            kk = np.polyfit(px, py, 1)
            pred, slope = np.polyval(kk, x), abs(kk[0])
        elif hist:
            pred, slope = hist[-1][1], 0.0
        else:
            pred, slope = None, 0.0
        if pred is None:
            y = max(cl)                    # start at the lowest blue cluster at the right edge (the trace)
        else:
            y = min(cl, key=lambda c: abs(c - pred))
            if abs(y - pred) > max_jump + 2.0 * slope:
                continue
        hist.append((x, y))
        path[x] = y
    return path


SPUR = (3.15e4, 3.9e4)       # Hz: narrow spur in the trace, excluded and bridged (the analyzer suppresses spurs too)


def track_lr(mask, x_start, x_stop, max_jump=6.0):
    """Left-to-right pass seeded with the topmost trace cluster in the first column."""
    hist, path = [], {}
    for x in range(x_start, x_stop):
        cl = clusters(mask[:, x])
        if not cl:
            continue
        if hist:
            pred = hist[-1][1] if len(hist) < 10 else np.polyval(
                np.polyfit([a for a, _ in hist[-20:]], [b for _, b in hist[-20:]], 1), x)
            y = min(cl, key=lambda c: abs(c - pred))
            if abs(y - pred) > max_jump:
                continue
        else:
            y = min(cl)
        hist.append((x, y))
        path[x] = y
    return path


def digitize(img):
    r0, r1, c0, c1 = PLOT
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV).astype(int)
    h, s, v = hsv[..., 0], hsv[..., 1], hsv[..., 2]
    m = (h >= 108) & (h <= 128) & (s > 50) & (v > 90)
    keep = np.zeros(m.shape, bool)
    keep[r0:r1, c0:c1] = True
    for a, bb, c, d in EXCLUDE:
        keep[a:bb, c:d] = False
    m &= keep
    right = track(m, c1 - 1, c0)                  # right edge -> left; stops at the ~35 kHz spur
    left = track_lr(m, c0 + 1, int(x_of(SPUR[0])))  # left edge -> right, up to the spur
    xs_spur = (x_of(SPUR[0]), x_of(SPUR[1]))
    path = {x: y for x, y in left.items() if x < xs_spur[0]}
    path.update({x: y for x, y in right.items() if x > xs_spur[1]})   # spur window bridged by interpolation
    xs = np.array(sorted(path), float)
    ys = np.array([path[int(x)] for x in xs])
    f, d = f_of(xs), db_of(ys)
    sel = (F_GRID >= f.min()) & (F_GRID <= f.max())
    pts = [[float(g), float(np.interp(np.log(g), np.log(f), d))] for g in F_GRID[sel]]
    return pts, (xs, ys)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pdf", required=True)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "gong_fig13b.json"))
    ap.add_argument("--overlay", default=None)
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as td:
        subprocess.run(["pdftoppm", "-r", "300", "-f", "10", "-l", "10", "-png", a.pdf, os.path.join(td, "p")], check=True)
        img = cv2.imread(glob.glob(os.path.join(td, "p-*.png"))[0])
    pts, (xs, ys) = digitize(img)
    res = {"source": "Gong et al., IEEE JSSC 58(5), 2023, Fig. 13(b): 4.2 K, PLL at ~10 GHz, measured at the /4 output (2.500018 GHz)",
           "units": "offset [Hz], SSB phase noise [dBc/Hz] at the 2.5 GHz test output; add 20log10(4)=12.04 dB for the 10 GHz VCO",
           "calibration": {"px_per_decade": DEC, "px_per_dB": float(-_k)},
           "instrument_reference": {"spot_dBc_Hz": {"1e4": -127.03, "1e5": -132.74, "1e6": -133.87, "1e7": -149.01},
                                    "integrated_10k_30M": {"int_noise_dBc": -67.80, "PM_urad": 575.92, "FM_Hz": 2308.0,
                                                           "jitter_fs": 36.663}},
           "trace": pts}
    json.dump(res, open(a.out, "w"), indent=1)
    print("wrote", a.out, len(pts), "points;", len(xs), "traced columns")
    if a.overlay:
        ov = img.copy()
        for x, y in zip(xs[::3], ys[::3]):
            cv2.circle(ov, (int(x), int(y)), 2, (0, 0, 255), -1)
        for dbv, yy in GRID_ROWS.items():
            cv2.line(ov, (219, int(round(y_of(dbv)))), (1240, int(round(y_of(dbv)))), (0, 200, 0), 1)
        cv2.imwrite(a.overlay, ov[900:1400, 200:1260])


if __name__ == "__main__":
    main()
