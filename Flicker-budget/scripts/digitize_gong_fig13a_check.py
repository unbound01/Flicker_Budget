#!/usr/bin/env python3
"""Second validation case: Gong et al. JSSC 2023 Fig. 13(a) (300 K, same analyzer screenshot layout as Fig. 13(b)),
traced with the Fig. 13(b) method but its own calibration. Its printed integrated values (10 kHz-30 MHz) are the check.

Usage: python scripts/digitize_gong_fig13a_check.py --pdf PATH/TO/Gong_JSSC2023.pdf [--out data/gong_fig13a.json] [--overlay out.png]
Output is derived from a copyrighted figure: keep it out of public repos.
"""
import argparse, glob, json, os, subprocess, sys, tempfile

import cv2
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import digitize_gong_fig13b as d13  # noqa: E402

d13.GRID_ROWS = {-135: 422.0, -140: 475.0, -145: 528.0, -150: 581.0, -160: 687.0, -165: 741.0}
d13._k, d13._b = np.polyfit(list(d13.GRID_ROWS), list(d13.GRID_ROWS.values()), 1)
d13.X_100K, d13.X_1M, d13.X_10M = 511.0, 805.0, 1099.0
d13.DEC = (d13.X_10M - d13.X_100K) / 2
d13.X_10K = d13.X_100K - d13.DEC
d13.PLOT = (312, 793, 219, 1239)
d13.EXCLUDE = []
d13.SPUR = (4.0e4, 5.2e4)
REF = {"spot_dBc_Hz": {"1e4": -125.77, "1e5": -130.36, "1e6": -128.79, "1e7": -142.48},
       "integrated_10k_30M": {"int_noise_dBc": -61.55, "PM_urad": 1180.0, "FM_Hz": 5949.0, "jitter_fs": 75.361}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", required=True)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "data", "gong_fig13a.json"))
    ap.add_argument("--overlay", default=None)
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as td:
        subprocess.run(["pdftoppm", "-r", "300", "-f", "10", "-l", "10", "-png", a.pdf, os.path.join(td, "p")], check=True)
        img = cv2.imread(glob.glob(os.path.join(td, "p-*.png"))[0])
    pts, (xs, ys) = d13.digitize(img)
    t = np.array(pts)
    for k, v in REF["spot_dBc_Hz"].items():
        print(f"  marker {float(k):8.0e} Hz: digitized {np.interp(np.log(float(k)), np.log(t[:, 0]), t[:, 1]):8.2f}  instrument {v}")
    json.dump({"source": "Gong et al., IEEE JSSC 58(5), 2023, Fig. 13(a): 300 K, /4 output (2.499917 GHz)",
               "instrument_reference": REF, "trace": pts}, open(a.out, "w"), indent=1)
    print("wrote", a.out)
    if a.overlay:
        ov = img.copy()
        for x, y in zip(xs[::3], ys[::3]):
            cv2.circle(ov, (int(x), int(y)), 2, (0, 0, 255), -1)
        cv2.imwrite(a.overlay, ov[250:890, 200:1260])


if __name__ == "__main__":
    main()
