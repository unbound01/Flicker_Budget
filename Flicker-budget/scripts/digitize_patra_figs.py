#!/usr/bin/env python3
"""Digitize Figs. 22, 24 and 25 of Patra et al., IEEE JSSC 53(1), 2018 (pp. 317-318).

Usage:
    python scripts/digitize_patra_figs.py --pdf PATH/TO/CryoCMOS_BPatra.pdf \
        --out data/patra_digitized.json [--overlays figures/]

You must supply your own copy of the paper; it is not part of this repository.

IMPORTANT: the pixel calibrations below were measured on the 13-page IEEE Xplore PDF
(PDF pages 9 and 10, rendered at 300 dpi -> 2550x3300 px). A different PDF layout needs
re-calibration (axis frames and gridlines).

Method: render the pages, locate axis frames/gridlines (constants below), extract curves
by colour mask (Fig. 22, tracked from the clean right edge leftwards so that label text and
arrows are ignored) and markers by morphological opening (Figs. 24, 25).
Estimated error: ~ +-0.3 dB per point (1 px = 0.15 dB in Fig. 22).
"""
import argparse
import glob
import json
import os
import subprocess
import tempfile

import cv2
import numpy as np

DPI = 300

# ---- Fig. 22 (PDF page 9; crop = right half, rows 8%..46%) -------------------------------
F22_CROP = (0.50, 0.08, 1.00, 0.46)       # fractions of (W, H, W, H)
F22_X0, F22_X1 = 170.5, 1019.5            # px of 1e5 Hz and 3e7 Hz (log axis)
F22_Y0, F22_Y1 = 452.5, 1049.5            # px of -80 dBc/Hz and -170 dBc/Hz
F22_MASKS = {                              # HSV colour masks of the three curves
    "4K":   lambda h, s, v: ((h < 8) | (h > 172)) & (s > 150) & (v > 150),
    "300K": lambda h, s, v: (h > 105) & (h < 135) & (s > 80) & (v > 90) & (v < 200),
    "165K": lambda h, s, v: (h > 45) & (h < 85) & (s > 120) & (v > 60),
}

# ---- Fig. 24 (PDF page 10): corner vs temperature, log-log ---------------------------------
F24_X1K, F24_PX_PER_DEC_X = 321.0, 355.0   # px of 1 K ; px per decade in T
F24_Y1E7, F24_PX_PER_DEC_Y = 892.5, 266.5  # px of 1e7 Hz ; px per decade in corner
F24_REGION = (897, 1420, 326, 1240)        # rows, cols

# ---- Fig. 25 (PDF page 10): PN at 4 K across tuning range, log x-axis ----------------------
F25_GRID_PX = [1444, 1577, 1705, 1829.5, 1950, 2067, 2181, 2291.0]
F25_GRID_GHZ = [5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2]
F25_Y_TOP_PX, F25_PX_PER_10DB = 254.0, 71.0   # px of -110 dBc/Hz ; px per 10 dB
F25_REGION = (256, 679, 1380, 2360)
F25_MASKS = {
    "1MHz":  lambda h, s: (h >= 98) & (h <= 110) & (s > 150),
    "10MHz": lambda h, s: (h >= 4) & (h <= 11) & (s > 150),
    "30MHz": lambda h, s: (h >= 17) & (h <= 25) & (s > 150),
}


def render(pdf, page, outdir):
    prefix = os.path.join(outdir, f"p{page}")
    subprocess.run(["pdftoppm", "-r", str(DPI), "-f", str(page), "-l", str(page), "-png", pdf, prefix],
                   check=True)
    return cv2.imread(glob.glob(prefix + "-*.png")[0])


def _clusters(col):
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


def _track(mask, x0, x1):
    """Follow a curve from the right edge leftwards, rejecting jumps (labels, arrows, spurs)."""
    hist, path = [], {}
    for x in range(int(x1) - 4, int(x0) + 4, -1):
        cl = _clusters(mask[:, x])
        if len(hist) >= 8:
            px = np.array([a for a, _ in hist[-30:]], float)
            py = np.array([b for _, b in hist[-30:]], float)
            pred = np.polyval(np.polyfit(px, py, 1), x)
        elif hist:
            pred = hist[-1][1]
        else:
            pred = None
        if not cl:
            continue
        if pred is None:
            y = cl[0]
        else:
            y = min(cl, key=lambda c: abs(c - pred))
            if abs(y - pred) > 7:
                continue
        hist.append((x, y))
        path[x] = y
    return path


def digitize_fig22(page9, overlay=None):
    H, W = page9.shape[:2]
    fx0, fy0, fx1, fy1 = F22_CROP
    im = page9[int(H * fy0):int(H * fy1), int(W * fx0):int(W * fx1)]
    h_, w_ = im.shape[:2]
    dec = (F22_X1 - F22_X0) / np.log10(3e7 / 1e5)
    pdb = (F22_Y1 - F22_Y0) / 90.0
    inv_x = lambda px: 1e5 * 10 ** ((px - F22_X0) / dec)
    inv_y = lambda py: -80 - (py - F22_Y0) / pdb
    hsv = cv2.cvtColor(im, cv2.COLOR_BGR2HSV)
    h, s, v = [hsv[..., i].astype(int) for i in range(3)]
    inside = np.zeros((h_, w_), bool)
    inside[int(F22_Y0) + 3:int(F22_Y1) - 2, int(F22_X0) + 3:int(F22_X1) - 2] = True
    grid = np.logspace(np.log10(1.5e5), np.log10(2.9e7), 40)
    out = {}
    for name, fn in F22_MASKS.items():
        path = _track(fn(h, s, v) & inside, F22_X0, F22_X1)
        px = np.array(sorted(path))
        f = inv_x(px)
        d = inv_y(np.array([path[p] for p in px]))
        out[name] = [[float(g), float(np.interp(np.log10(g), np.log10(f), d))]
                     for g in grid if f.min() <= g <= f.max()]
        if overlay is not None:
            for g, dd in out[name]:
                cv2.circle(overlay, (int(F22_X0 + np.log10(g / 1e5) * dec) + int(W * fx0),
                                     int(F22_Y0 + (-80 - dd) * pdb) + int(H * fy0)), 4, (0, 0, 255), 1)
    return out


def _markers(mask, k):
    op = cv2.morphologyEx(mask.astype(np.uint8), cv2.MORPH_OPEN,
                          cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k)))
    n, _, st, cen = cv2.connectedComponentsWithStats(op)
    return sorted((cen[i][0], cen[i][1]) for i in range(1, n) if st[i][4] > 20)


def digitize_fig24(page10, overlay=None):
    hsv = cv2.cvtColor(page10, cv2.COLOR_BGR2HSV)
    h, s, v = [hsv[..., i].astype(int) for i in range(3)]
    r0, r1, c0, c1 = F24_REGION
    reg = np.zeros(h.shape, bool)
    reg[r0:r1, c0:c1] = True
    blue = (h > 100) & (h < 135) & (s > 120) & (v > 70) & reg
    pts = _markers(blue, 15)
    out = [{"T_K": float(10 ** ((x - F24_X1K) / F24_PX_PER_DEC_X)),
            "corner_Hz": float(10 ** (7 - (y - F24_Y1E7) / F24_PX_PER_DEC_Y))} for x, y in pts]
    if overlay is not None:
        for x, y in pts:
            cv2.circle(overlay, (int(x), int(y)), 16, (0, 0, 255), 3)
    return out


def digitize_fig25(page10, overlay=None):
    hsv = cv2.cvtColor(page10, cv2.COLOR_BGR2HSV)
    h, s = hsv[..., 0].astype(int), hsv[..., 1].astype(int)
    r0, r1, c0, c1 = F25_REGION
    reg = np.zeros(h.shape, bool)
    reg[r0:r1, c0:c1] = True
    b, a = np.polyfit(np.log10(F25_GRID_GHZ), F25_GRID_PX, 1)     # px = a + b*log10(GHz)
    out = {}
    for name, fn in F25_MASKS.items():
        pts = _markers(fn(h, s) & reg, 9)
        rec = [{"f_GHz": float(10 ** ((x - a) / b)),
                "PN": float(-110 - (y - F25_Y_TOP_PX) / F25_PX_PER_10DB * 10)} for x, y in pts]
        out[name] = [r for r in rec if abs(r["f_GHz"] - 7.011) > 0.02]    # drop legend markers
        if overlay is not None:
            for r in out[name]:
                cv2.circle(overlay, (int(a + b * np.log10(r["f_GHz"])),
                                     int(F25_Y_TOP_PX + (-110 - r["PN"]) / 10 * F25_PX_PER_10DB)),
                           14, (0, 160, 0), 3)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pdf", required=True, help="your copy of the Patra et al. PDF")
    ap.add_argument("--out", default="data/patra_digitized.json")
    ap.add_argument("--overlays", default=None,
                    help="optional directory for overlay images (they contain the paper's figures: do not publish)")
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as td:
        p9, p10 = render(a.pdf, 9, td), render(a.pdf, 10, td)
    ov9 = p9.copy() if a.overlays else None
    ov10 = p10.copy() if a.overlays else None
    result = {
        "source": "Patra et al., JSSC 53(1), 2018, Figs 22, 24, 25 (300 dpi render, colour-masked / marker-detected)",
        "digitization_error_estimate": "~ +-0.3 dB per point; dominant uncertainty is inter-figure spread (~0.1-1.8 dB)",
        "fig22_PN_dBc_per_Hz_vs_offset_Hz": digitize_fig22(p9, ov9),
        "fig24_corner_Hz_vs_T_K": digitize_fig24(p10, ov10),
        "fig25_PN_dBc_per_Hz_vs_GHz": digitize_fig25(p10, ov10),
    }
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(result, f, indent=1)
    print("wrote", a.out)
    if a.overlays:
        os.makedirs(a.overlays, exist_ok=True)
        cv2.imwrite(os.path.join(a.overlays, "overlay_fig22.png"), ov9)
        cv2.imwrite(os.path.join(a.overlays, "overlay_fig24_25.png"), ov10)


if __name__ == "__main__":
    main()
