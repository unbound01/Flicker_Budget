# Phase 2: how well is the fitted flicker corner constrained?

Script: `scripts/fit_identifiability.py` (needs the locally generated data; see `data/README.md`).
Model: L(f) = α/f²·(1 + f_c/f), fitted in dB to the 4 K curve of Patra Fig. 22. Primary range 200 kHz-10 MHz.
The figure `figures/fit_sensitivity.png` is generated locally and git-ignored.

## Result in one paragraph

The **statistical** uncertainty of α and f_c is small (10-14%), but their **systematic** spread over reasonable fit ranges is
large (α ×3.3, f_c ×5.4), and the residuals show a smooth ±1 dB curvature that the model cannot follow. So α and f_c are
**poorly constrained individually (strongly range-dependent)** by this curve. Their product K = α·f_c is (±25% across ranges), and the quantity the
PLL budget cares about, FN with a 1 MHz PLL, is robust: 2.11-2.56 kHz in every fit, never below 1.9 kHz.

## Parameter-uncertainty table (4 K, primary fit unless stated)

| Quantity | Primary fit | Linearised 1σ (σ = 0.3 dB / 0.63 dB) | Profile likelihood 68% (95%) | MC, independent 0.3 dB (16-84%) | MC, all errors (16-84%) | 16 fit ranges (min-max) |
|---|---|---|---|---|---|---|
| α (Hz²/Hz) | 0.137 | ×/÷ 1.05 / 1.10 | n/a | 0.130-0.144 | 0.101-0.181 | 0.072-0.236 |
| f_c (MHz) | 4.12 | ×/÷ 1.06 / 1.14 | 3.61-4.75 (3.18-5.44) | 3.86-4.43 | 2.81-6.06 | 1.96-10.50 |
| K = α·f_c (10⁵ Hz²) | 5.66 | ×/÷ 1.02 / 1.05 | n/a | 5.52-5.80 | 5.01-6.24 | 4.62-7.60 |
| Free-running FN (kHz) | 3.28 | n/a | n/a | 3.26-3.30 | 3.21-3.36 | n/a |
| FN, 1 MHz PLL, fixed floor (kHz) | 2.27 | n/a | n/a | 2.25-2.29 | 2.16-2.41 | 2.11-2.56 |
| PLL bandwidth for 1.9 kHz (MHz) | 2.90 | n/a | n/a | 2.75-3.06 | 2.11-3.99 | n/a |

"All errors" = independent 0.2 dB + common offset 0.2 dB + a tilt of ±1 dB at 29 MHz (the Fig. 22 / Fig. 25 disagreement), 1500 refits per scenario.
Probability of FN (1 MHz PLL) below 1.9 kHz: 0 in all 6000 samples.

## Findings

1. **The model does not fully fit the data.** Residual rms 0.61 dB against ~0.3 dB digitization error (χ²_red = 4.5 at 0.3 dB). Residuals versus
   offset are U-shaped in every fit (about +1 dB at the range edges, −0.6 dB in the middle, growing to +2.5 dB at 29 MHz outside the 10 MHz range), so the
   parameters drift with the range: f_c grows with both the lower and the upper fit edge.
2. **α and f_c trade off.** corr(ln α, ln f_c) = −0.96; the ridge slope is −1.33 (−1 would be constant K). Data pin down K (±25%) far better than α or f_c.
3. **Why.** The flicker term exceeds the white term below f_c, which is 77% of the fitted log-frequency span; the white-dominated tail
   (about 4-10 MHz) is short and sits where the curve deviates most from the model, so α is weakly constrained and f_c = K/α inherits that.
4. **Comparison with Fig. 24 (5.68 MHz at 4 K).** For the primary range the 95% profile interval (3.2-5.4 MHz) excludes it (Δχ² = 4.9); fixing f_c = 5.68 MHz
   raises the rms only from 0.61 to 0.67 dB and K by 10%. Across the range grid the fitted f_c spans 2.0-10.5 MHz and includes 5.7. The comparison is **not decisive**.
5. **PLL conclusion (brick-wall model).** α and f_c move along the K ridge in compensating directions, so FN with a 1 MHz PLL changes little (2.11-2.56 kHz over windows, 2.16-2.41 kHz over
   digitization error). None reaches 1.9 kHz, so the conditional conclusion survives under the brick-wall model; the generic-loop results are in `pll_model.md`.

## Justified fit range

Keep **200 kHz-10 MHz** as the reference range: the lower edge clears the supply spurs near 130 kHz, and the upper edge is Patra's f_b and stays below the
10-30 MHz region where Figs. 22 and 25 disagree by up to 1.9 dB. Do **not** adjust the range to match Fig. 24: 300 kHz gives 5.1 MHz and 500 kHz gives 6.5 MHz, so that would be circular.
Report K and the PLL-predicted FN as the robust results, and quote the grid spread as the systematic uncertainty on α and f_c.

## How to read these intervals

- The linearised, profile-likelihood and Monte Carlo intervals are **conditional on the two-term model** and treat the 29 fitted points as independent. They are not:
  the points are resampled from one traced curve, so neighbours share pixels and errors are correlated (the Monte Carlo includes common-offset and tilt components, which only partly
  covers this). Treat them as sensitivity measures, not calibrated 95% confidence intervals.
- The 16-window grid is a hand-chosen **systematic sensitivity study**, not a probability distribution.
- "Poorly constrained" means the fitted values move a lot with the window and the model form; it does not prove the parameters are mathematically non-identifiable.
- The 165 K and 300 K fits (f_c ≈ 0.17 and 0.05 MHz, at or below the data range) are weakly constrained, so any K-versus-temperature ratio is exploratory.

## Not done

- A model with a free slope or a physical oscillator spectrum, which the U-shaped residuals call for.
- Data below 100 kHz, or the authors' raw curves.
- ~~A PLL transfer-function model~~: done as a generic linear model in Phase 3 (`docs/pll_model.md`).
