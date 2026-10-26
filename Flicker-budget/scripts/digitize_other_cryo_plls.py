#!/usr/bin/env python3
"""Digitize the measured 4 K phase noise of two independent cryo-CMOS PLLs (not from the Delft group):

  Xiao et al., "A 19.4-fs_RMS Jitter 0.1-to-44GHz Cryo-CMOS Fractional-N CP-PLL ...", IEEE RFIC 2025,
      DOI 10.1109/RFIC61188.2025.11082927, Fig. 7 (top right): integer mode, 8.6 GHz, analyzer trace T1 = 4 K.
      Targets printed in the screenshot: spot noise 1 kHz-100 MHz and integrated 12 kHz-20 MHz
      (-62.60 dBc, PM 1.05 mrad, FM 3.333 kHz, 19.485 fs).
  Peng et al., "A 94.4fs_rms Cryogenic Sub-Sampling PLL in 4K ...", IEEE RFIT 2024,
      DOI 10.1109/RFIT60557.2024.10812563, Fig. 6 (a)-(d): integer-N, 5.7/6.0/6.3/6.9 GHz, blue = 4 K.
      Targets printed in each panel: PN at 100 kHz and jitter 1 kHz-100 MHz.

Both figures are low-resolution bitmaps in the PDFs (about 350 px wide), so expect ~0.5 dB per point.
Usage:
  python scripts/digitize_other_cryo_plls.py --xiao PATH/Xiao_RFIC2025.pdf --peng PATH/Peng_RFIT2024.pdf
      [--out-xiao data/xiao_rfic2025_fig7.json] [--out-peng data/peng_rfit2024_fig6.json] [--overlay-dir DIR]
Outputs are derived from copyrighted figures: keep them out of public repos.
Calibration is for the IEEE Xplore PDFs, page 3 rendered at 600 dpi (5100 x 6600 px).
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
from digitize_gong_thesis_fig721 import track_robust  # noqa: E402


def render(pdf, page, td, tag):
    subprocess.run(["pdftoppm", "-r", "600", "-f", str(page), "-l", str(page), "-png", pdf, os.path.join(td, tag)], check=True)
    return cv2.imread(glob.glob(os.path.join(td, tag + "-*.png"))[0])


def smooth(ys, k=3):
    return np.array([np.median(ys[max(0, i - k):i + k + 1]) for i in range(len(ys))])


def xiao(img):
    crop = img[396:1419, 3748:5023]                    # Fig. 7 top-right screenshot
    hsv = cv2.cvtColor(crop, cv2.COLOR_BGR2HSV).astype(int)
    h, s, v = hsv[..., 0], hsv[..., 1], hsv[..., 2]
    m = (h >= 108) & (h <= 128) & (s > 45) & (v < 245)  # navy 4 K trace (red = 300 K, magenta = estimate)
    keep = np.zeros(m.shape, bool)
    keep[160:775, 2:952] = True
    for r0, r1, c0, c1 in ((300, 410, 310, 640), (505, 610, 310, 640), (140, 330, 735, 960), (685, 750, 640, 960)):
        keep[r0:r1, c0:c1] = False                       # overlaid labels, spot-noise list, "PN Estimation"
    path = track_robust(m & keep, 3, 951, max_jump=8, reacquire=10)
    xs = np.array(sorted(path), float)
    ys = smooth(np.array([path[int(x)] for x in xs]))
    a, b = np.polyfit([182, 375, 567, 763.], [4, 5, 6, 7.], 1)          # decade gridlines -> log10 f
    k, c = np.polyfit([206, 264, 318, 376, 435, 492, 549, 606, 660, 717.],
                      [-80, -90, -100, -110, -120, -130, -140, -150, -160, -170.], 1)
    f, d = 10 ** (a * xs + b), k * ys + c
    grid = np.logspace(np.log10(max(f[0], 1.2e3)), np.log10(min(f[-1], 9.9e7)), 200)
    pts = [[float(g), float(np.interp(np.log(g), np.log(f), d))] for g in grid]
    out = {"source": "Xiao et al., IEEE RFIC 2025, Fig. 7, integer mode, 8.6 GHz, 4 K (trace T1)", "f0": 8.599999936e9,
           "instrument_reference": {"spot": {"1e3": -106.69, "1e4": -116.14, "1e5": -123.06, "1e6": -129.39,
                                             "1e7": -145.26, "1e8": -156.27},
                                    "spot_note": "read from a blurred bitmap; the 1 MHz digits are uncertain",
                                    "int_12k_20M": {"int_noise_dBc": -62.60, "PM_urad": 1050.0, "FM_Hz": 3333.0,
                                                    "jitter_fs": 19.485}},
           "trace": pts}
    return out, crop, (xs, ys)


PENG_PANELS = {"5.7GHz": (212.5, 1028, 107.5, 582.5, -114.2194, 95.8, 5.7e9),
               "6.0GHz": (1215.5, 2031, 107.5, 582.5, -113.6893, 94.4, 6.0e9),
               "6.3GHz": (212.5, 1028, 754.5, 1230.5, -113.8069, 97.7, 6.3e9),
               "6.9GHz": (1215.5, 2031, 754.5, 1230.5, -113.8512, 97.2, 6.9e9)}


def peng(img):
    crop = img[396:2640, 306:2448]                     # Fig. 6, four MATLAB panels, frames 1 kHz-100 MHz x -160..-90 dBc/Hz
    hsv = cv2.cvtColor(crop, cv2.COLOR_BGR2HSV).astype(int)
    h, s, v = hsv[..., 0], hsv[..., 1], hsv[..., 2]
    blue = (h >= 108) & (h <= 130) & (s > 110) & (v > 70)
    out = {"source": "Peng et al., IEEE RFIT 2024, Fig. 6, integer-N, 4 K curves", "panels": {}}
    marks = []
    for name, (xl, xr, yt, yb, pn100, jit, f0) in PENG_PANELS.items():
        X = lambda f: xl + (np.log10(f) - 3) * (xr - xl) / 5
        Y = lambda dd: yt + (-90 - dd) * (yb - yt) / 70
        keep = np.zeros(blue.shape, bool)
        keep[int(yt) + 3:int(yb) - 2, int(xl) + 3:int(xr) - 2] = True
        keep[int(Y(-133.5)):int(Y(-156)), int(xl):int(X(4.5e6))] = False      # blue "4K ..." text block
        keep[int(yt):int(Y(-102)), int(X(1.2e7)):int(xr)] = False             # legend
        path = track_robust(blue & keep, int(xl) + 4, int(xr) - 3, max_jump=8, reacquire=10)
        xs = np.array(sorted(path), float)
        ys = smooth(np.array([path[int(x)] for x in xs]), 2)
        f = 10 ** (3 + (xs - xl) * 5 / (xr - xl))
        d = -90 - (ys - yt) * 70 / (yb - yt)
        grid = np.logspace(np.log10(f[0]), np.log10(f[-1]), 200)
        out["panels"][name] = {"f0": f0, "instrument_reference": {"pn100k": pn100, "jitter_1k_100M_fs": jit},
                               "trace": [[float(g), float(np.interp(np.log(g), np.log(f), d))] for g in grid]}
        marks.append((xs, ys))
    return out, crop, marks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xiao", required=True)
    ap.add_argument("--peng", required=True)
    ap.add_argument("--out-xiao", default=os.path.join(HERE, "..", "data", "xiao_rfic2025_fig7.json"))
    ap.add_argument("--out-peng", default=os.path.join(HERE, "..", "data", "peng_rfit2024_fig6.json"))
    ap.add_argument("--overlay-dir", default=None)
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as td:
        ix, ip = render(a.xiao, 3, td, "x"), render(a.peng, 3, td, "p")
    ox, cx, mx = xiao(ix)
    op, cp, mp = peng(ip)
    json.dump(ox, open(a.out_xiao, "w"), indent=1)
    json.dump(op, open(a.out_peng, "w"), indent=1)
    print("wrote", a.out_xiao, "and", a.out_peng)
    if a.overlay_dir:
        os.makedirs(a.overlay_dir, exist_ok=True)
        for crop, sets, name in ((cx, [mx], "overlay_xiao.png"), (cp, mp, "overlay_peng.png")):
            ov = crop.copy()
            for xs, ys in sets:
                for x, y in zip(xs[::3], ys[::3]):
                    cv2.circle(ov, (int(x), int(y)), 2, (0, 200, 0), -1)
            cv2.imwrite(os.path.join(a.overlay_dir, name), ov)


if __name__ == "__main__":
    main()
