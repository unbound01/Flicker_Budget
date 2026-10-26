# Phase 1: baseline and source validation

Purpose: show that the implementation reproduces Patra et al.'s original (white-FM-only) calculation
and that the digitized data are trustworthy before the flicker term is used for any stronger claim.
Numbers below come from `src/budget.py` and from data generated locally with `scripts/digitize_patra_figs.py`.

## 1. Source parameters (Patra et al., JSSC 53(1), 2018, pp. 311, 315, 317-318)

| Item | Value in the source |
|---|---|
| Target | LO frequency noise ≤ 1.9 kHz rms, quoted for a 6 GHz carrier |
| Band | 8.3 kHz (1/T₂*, T₂* = 120 µs) to 10 MHz (qubit operation speed) |
| Definition | L(Δf): single-sideband phase noise, dBc/Hz; FN² = 2∫Δf² L dΔf (Eq. 8) |
| Oscillator | L_osc = α/Δf², α in Hz²/Hz (beyond f_PLL) |
| PLL model | flat in-band floor equal to L_osc(f_PLL), f_PLL = 300 kHz "typical" (Eq. 7); 1/f³ corner assumed below f_PLL (footnote 3) |
| Quoted result | α < 0.2 Hz²/Hz; in-band −115 dBc/Hz; −147 dBc/Hz at 10 MHz (Eq. 9) |
| Measurements | Figs. 22, 24, 25 at 6.3 GHz; free-running FN 3.4 kHz (8.3 kHz to 10 MHz); 4 K corner ≈ 5.7 MHz (Fig. 24) |

## 2. Baseline reproduction (f_c = 0)

| Quantity | Patra | This code | Verdict |
|---|---|---|---|
| α bound | < 0.2 | 0.1842 exact; α = 0.2 gives 1.980 kHz (+4.2% over 1.9) | quoted bound is loose by 8% |
| PN at 10 MHz | −147 | −147.0 (α = 0.2), −147.35 (α_max) | reproduced |
| In-band PN at 300 kHz | −115 | −116.5 (α = 0.2), −116.9 (α_max) | **not reproduced, 1.5-1.9 dB; unexplained** |
| Closed form vs numerical integral | n/a | worst relative difference 3×10⁻¹⁶ | unit test of the algebra, not a result |
| Free-running FN, pure 1/f² | 3.4 kHz measured | α_eff = 0.58 would be needed (3× the budget) | motivates the 1/f³ term |

The in-band gap changes α_max by under 0.6%, so it does not affect any later result. Gate: **passed**, with the two deviations above recorded.

## 3. Digitization checks

| Check | Result |
|---|---|
| Fig. 22 frequency axis: gridlines at 1e6 and 2e7 Hz (not used for calibration) | −0.7 px and +0.4 px, i.e. 0.5% and 0.2% in frequency |
| Fig. 22 level axis: eight dotted gridlines, −90 to −160 dBc/Hz | mean error 0.09 dB, worst 0.18 dB |
| Fig. 24 temperature axis | 11 markers fall at 4.0, 20.1, 60.1, 70.3, 90.1, 125.6, 170.4, 200.3, 230.5, 250.6, 300.5 K |
| Overlay of extracted points on the original figures | points sit on the curves and markers (visual) |
| Fig. 22 (4 K) vs Fig. 25 at 6.3 GHz | 1 MHz: −0.66 dB; 10 MHz: −0.06 dB; 29 MHz: +1.86 dB |
| Re-extraction from the PDF | identical output (reproducibility, not accuracy) |

Estimated uncertainty: about ±0.3 dB per point (axis ≤ 0.2 dB plus curve-centre error; my estimate). The larger uncertainty is the **1.9 dB spread between figures at 30 MHz**, probably a measurement-floor or setup difference. Fig. 25's 30 MHz points were not used in any fit.

## 4. Fit range

Primary range 200 kHz to 10 MHz: the lower edge clears the supply spurs near 130 kHz in the 4 K curve, and the upper edge is Patra's f_b. It was chosen for those reasons, not optimized. Sensitivity (4 K):

| Fit | α (Hz²/Hz) | f_c (MHz) | rms (dB) | Free-running FN (kHz) | FN, 1 MHz PLL, fixed / continuous (kHz) |
|---|---|---|---|---|---|
| 200 kHz-10 MHz | 0.137 | 4.12 | 0.61 | 3.28 | 2.27 / 2.35 |
| 300 kHz-10 MHz | 0.122 | 5.12 | 0.50 | 3.36 | 2.27 / 2.36 |
| 150 kHz-29 MHz | 0.092 | 6.48 | 0.82 | 3.21 | 2.11 / 2.20 |
| 200 kHz-10 MHz, +0.5 dB | 0.154 | 4.12 | 0.61 | 3.48 | 2.41 / 2.49 |
| 200 kHz-10 MHz, −0.5 dB | 0.122 | 4.12 | 0.61 | 3.10 | 2.15 / 2.22 |

The fitted corner depends on the range; a full 16-range grid (`docs/identifiability.md`) widens the span to 2.0-10.5 MHz, which includes the Fig. 24 value (5.7 MHz), although the primary range's 95% profile interval (3.2-5.4 MHz) excludes it. The 1 MHz-PLL result stays above 1.9 kHz in every fit (2.11-2.56 kHz over the grid).

## 5. Not yet done

- A side-by-side read of Patra's Section V against §1 and §2 above by the author of this repo. Eqs. 7–9, footnotes 3 and 4 and Fig. 24 were checked against the rendered page on 16 Oct (Eqs. 8 and 9 carry a square root that the text extraction drops).
- A Fig. 23 cross-check (PN at 30 MHz versus temperature) against the 4 K curve. Partial check (16 Oct): at 4 K, Fig. 23 and Fig. 25 give about −162 dBc/Hz at 30 MHz, while the two-term fit gives −157.6 dBc/Hz, so the fit overestimates the PN there by about 4–5 dB. The fit should not be extrapolated beyond 10 MHz.
- Per-point validation against raw data, which requires the authors' measurements.
- ~~A PLL noise-transfer model~~: done as a generic linear model in Phase 3 (`docs/pll_model.md`); an architecture-specific model of Gong's PLL is Phase 4.
