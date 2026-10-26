# Candidate 8: tolerance of a PLL-based qubit LO to a 1/f³ corner

Week 3 report · 12 Oct 2026

## 1. Aim

Patra et al. (JSSC 2018) budget the LO frequency noise (FN) at ≤ 1.9 kHz rms over 8.3 kHz to 10 MHz, assuming the oscillator's 1/f³ corner lies below the PLL bandwidth (their footnote 3). We add a 1/f³ term, find the largest tolerable corner, and test the model on Patra's measured 4 K data.

## 2. Model and assumptions

- Oscillator: L(Δf) = α/Δf² · (1 + f_c/Δf), single-sideband, α in Hz²/Hz.
- PLL: brick wall at f_PLL = 300 kHz with a flat in-band floor, in two variants:
  - **fixed floor:** α/f_PLL²
  - **continuous floor:** α(1 + f_c/f_PLL)/f_PLL²
- FN² = 2∫Δf² L dΔf from f_a = 8.3 kHz (1/T₂*, T₂* = 120 µs) to f_b = 10 MHz.
- Closed form, with ξ = 1 − (f_a/f_PLL)³ ≈ 1 and K = α·f_c: FN² = 2α[(f_b − f_PLL) + ξ f_PLL/3] + 2K[ln(f_b/f_PLL) + ξ/3 (continuous floor only)]. The budget therefore depends only on α and K, and f_c,max = K_max(α)/α.
- Origin of 1.9 kHz: the qubit's own frequency noise for T₂* = 120 µs (van Dijk et al., Phys. Rev. Applied 12, 044054 (2019), Table II, which lists "nuclear spin noise 1.9 kHz rms, T₂* = 120 µs"; T₂* = 1/(√2 π σ_f) = 118 µs). It is a "do not degrade T₂*" heuristic, not an oscillator limit derived from fidelity.

## 3. Units and conventions (checked)

- All frequencies in Hz, not rad/s (σ_ω = 2π σ_f).
- L is single-sideband. α = L·Δf²; the one-sided FN spectrum is S_f = 2α, which gives the factor 2 in FN² (Patra Eq. 8). Ball et al. (Eq. 5, SI Eqs. 80–83) use the same single-sideband to one-sided relation, and their X-gate filter function (SI Eq. 112) gives 1−F = 9.06×10⁻⁷ for white FM at α = 0.184 and f_R = 1 MHz, against 9.08×10⁻⁷ from van Dijk and Gong.
- van Dijk Table II (published version) reproduced: L(1 MHz) = −106 dBc/Hz with ENBW 2.5 MHz gives 11.2 kHz (table: 11 kHz). σ_f = 1.9 kHz at f_R = 1 MHz gives 1−F = 3.61×10⁻⁶ (table: 3.6×10⁻⁶).
- For a π pulse, 1−F = (σ_f/f_R)² and ENBW = (π²/4)·f_R ≈ 2.47 f_R, so f_b = 10 MHz corresponds to f_R ≈ 4 MHz.
- α refers to the oscillator carrier and scales by N² if the carrier is multiplied. Patra's budget is for 6 GHz; Fig. 22 is at 6.3 GHz.

## 4. Results

**Baseline (f_PLL = 300 kHz).** α_max = 0.1842 Hz²/Hz (Patra: < 0.2; α = 0.2 gives 1.98 kHz). In-band PN −116.9 dBc/Hz (Patra −115; the 2 dB gap is unresolved). Patra's footnote 4 says −115 dBc/Hz is "easily achieved" by an integer-N sub-sampling PLL with a 50 MHz reference and a −155 dBc/Hz floor; referred to 6 GHz with 20·log₁₀(120) that floor is −113.4 dBc/Hz, above both −115 and −116.9. This assumes the −155 dBc/Hz is at the 50 MHz reference plane, which the footnote does not state. PN at 10 MHz −147.3 dBc/Hz (Patra −147).

**Corner tolerance f_c,max (kHz).**

| α (Hz²/Hz) | fixed floor | continuous floor |
|---|---|---|
| 0.17 | 233 | 213 |
| 0.15 | 637 | 582 |
| 0.10 | 2353 | 2149 |

At α_max with f_c = 500 kHz, FN = 2.063 kHz (white out-of-band 84%, flicker 15%, in-band 0.9%).

**Upper band limit (f_PLL = 300 kHz).**

| f_b (MHz) | 1 | 2 | 5 | 10 |
|---|---|---|---|---|
| α_max (Hz²/Hz) | 2.26 | 1.00 | 0.38 | 0.18 |
| equivalent f_R (MHz) | 0.41 | 0.81 | 2.03 | 4.05 |

**Gong 2022 oscillator (4.2 K, β = 0.041, f_C = 1.3 MHz).** FN = 1.09 kHz (fixed) or 1.10 kHz (continuous); f_c,max = 9.8 or 8.9 MHz, about 7× margin.

**Fit to Patra's 4 K data (Fig. 22, 6.3 GHz, 200 kHz to 10 MHz).**

- α = 0.137 Hz²/Hz, f_c = 4.1 MHz, rms residual 0.61 dB. Over a grid of 16 fit ranges α = 0.07–0.24 and f_c = 2.0–10.5 MHz, so the two are poorly constrained (strongly range-dependent); K = α·f_c varies only 0.46–0.76×10⁶ Hz² (docs/identifiability.md). Fig. 24 gives 5.7 MHz at 4 K; the primary range's 95% profile interval (3.2–5.4 MHz) excludes it, but the range-to-range spread is larger, so the comparison is not decisive.
- Free-running FN: 3.1–3.5 kHz from the fit, 3.26 kHz from direct integration of the digitized curve; Patra reports 3.4 kHz. The single (α, f_c) pair reproduces the 3.4 kHz and follows the PN curve to 0.61 dB rms, with a systematic curvature (see Limitations).
- Through a PLL: 2.11–2.56 kHz at 1 MHz in all 16 fit ranges (none reaches 1.9 kHz) and 2.4–2.8 kHz at 300 kHz (primary fits). For the primary fit, 1.9 kHz needs 2.9–3.6 MHz (fixed and continuous floors), or 2.1–4.0 MHz (68%) with digitization error. At α = 0.137, f_c,max is 0.96 MHz (300 kHz PLL) and 1.66 MHz (1 MHz PLL), against a fitted 4.1 MHz.
- Patra's conclusion that a PLL of ≤ 1 MHz can easily reach 1.9 kHz is not reproduced for this oscillator under the extended spectrum and the tested PLL/floor models. A generic linear phase-domain PLL model (docs/pll_model.md) shows the brick wall is an optimistic bound: with Patra's own white-FM floor α/f² the generic first-order, type-II and second-order loops need 4.6–6.0 MHz (4.4–5.8 MHz from the direct digitized spectrum) instead of 2.9 MHz, and 5.9–7.4 MHz with a flicker-aware floor instead of 3.6 MHz (3.1–8.4 MHz over the 16 fit windows). Under an illustrative −122 dBc/Hz floor scenario (Gong's measurement scaled to 6.3 GHz) no bandwidth reaches 1.9 kHz (minimum 2.8 kHz). At 1 MHz these models give 2.7 kHz (white floor) to 4.1 kHz, against 3.28 kHz free-running. These are scenario results for generic loops, not a proof that no PLL can meet the budget. The generic family closest to Gong's DAPLL (type-II, ζ = 1.6) needs 5.09 MHz (white-FM floor) or 6.43 MHz (flicker-aware floor), and at Gong's own 4 MHz gives 2.06 kHz with the white-FM floor.
- Across the tuning range (Fig. 25, from the 1 and 10 MHz points), FN with a 1 MHz PLL is 1.0–1.35 kHz at 5.7–6.05 GHz, 2.1–2.5 kHz at 6.2–6.35 GHz and 3.3–4.9 kHz at 6.55–7.3 GHz.

**K = α·f_c versus temperature.** K is 4.6–7.6×10⁵ Hz² at 4 K (16 fit ranges) and 2–7×10⁴ Hz² at 300 K, a rise of roughly 6–50× on cooling. Exploratory only: the 165 K and 300 K fits are weakly constrained, so this ratio is not a measured result. Gong's oscillator showed roughly constant flicker. The 300 K corner is poorly constrained because the data start at 150 kHz. The literature is split: Lei (citing Mehrpoo 2019) and Patra report higher flicker at cryogenic temperature, Gong 2022 reports no temperature dependence.

## 5. Verification

- Closed form versus numerical integration (scipy quad): worst relative difference 3×10⁻¹⁶ over α, f_c, f_PLL and both floors. This is a **unit test of the algebra, not a result**.
- Digitization: points overlaid on the original figures; Fig. 24 temperatures fall on round values. Error about ±0.3 dB per point. Fig. 22 and Fig. 25 differ at 6.3 GHz by 0.7 dB (1 MHz), 0.1 dB (10 MHz) and 1.8 dB (29 MHz).
- Independent route: Fig. 25 points at 6.2–6.35 GHz give a free-running FN of 3.2–4.0 kHz.

## 6. Limitations

- **Brick-wall PLL.** Real loops have a transition region and peaking; the linear model of Phase 3 uses generic first-, second-order and type-II loops with one aggregated flat input floor sent through a single transfer function, no cross-spectra, and no circuit-specific detail (extra filtering, shaped floors, injection locking, separate injection points). The in-band floor is not flat: reference and input-buffer flicker dominate below about 100 kHz (Gong 2023).
- **Discontinuity at f_PLL.** The fixed floor jumps by (1 + f_c/f_PLL) at f_PLL, which is 2.67 at f_c = 500 kHz and f_PLL = 300 kHz. This is a modelling choice, not physics; the continuous variant is the robustness check, and the two differ by about 10% in f_c,max.
- **Digitized data.** ±0.3 dB per point, with larger spread between figures. Fig. 22 starts at 100 kHz, so about 31% of the 4 K free-running FN² is extrapolated (the PLL results use only data above 300 kHz plus the floor assumption). Corners from Fig. 25 are ill-conditioned (3–18 MHz), so only its FN values are used. The fitted α and f_c depend strongly on the fit range (the residuals show a systematic ±1 dB curvature the model cannot follow); only K and the PLL-predicted FN are robust. The 165 K and 300 K fits are not well constrained.
- **Budget, not measurement.** 1.9 kHz is a heuristic, not a fidelity limit. For white FM, 99.999% fidelity at f_R = 1 MHz allows σ_f ≈ 3.2 kHz, so 1.9 kHz is conservative; with a flicker term the outcome depends on f_C (Gong 2022).
- **Lower cutoff.** f_a = 8.3 kHz versus an echo-period cutoff: a 10× change in f_a shifts the free-running fitted corner by about 50%. For an f⁻³ term, the Ramsey and primitive-gate filters (∝ ω² at low offset) give a logarithmic low-frequency divergence, while the spin-echo filter (∝ ω⁴) converges; the cutoff therefore matters for free evolution and driven gates, not for echo (Ball SI Eqs. 89, 90, 112).
- **Carrier.** Fits are at 6.3 GHz; the budget is quoted at 6 GHz.
- **Fit range.** The two-term fit is used only up to 10 MHz. At 30 MHz it gives −157.6 dBc/Hz against about −162 dBc/Hz in Patra's Figs. 23 and 25 (4 K), so it must not be extrapolated to wider bands such as Gong's 10 kHz–30 MHz jitter band.
- **Open items.** The 2 dB in-band gap remains open, and Patra's attribution of the lower band limit to the echo period is only partly supported: Ball's supplement gives no explicit lower bound, only the spin-echo filter, which scales as ω⁴ at low offset.

## 7. Prior art

Gong 2022 (IEEE TCAS-I 69(12)) gives β versus f_C under a fidelity criterion without a PLL, and notes that f_C can rise to about 10 MHz at 4.2 K. Gong 2023 (JSSC 58(5)) is the closest qualitative precedent: its specification uses the same brick-wall form as Patra (Eq. 3) and states that it neglects flicker (footnote 2); it notes that the VCO "would exhibit a very high flicker corner (e.g., a few MHz) at 4.2 K" that "could be partially mitigated by using a large DAPLL loop bandwidth" (p. 1370); and its measured 4.2 K profile is limited from 100 kHz to 4 MHz by filtered VCO flicker and INBUF shot noise (p. 1371). What is not found there is a quantitative tolerance: f_c,max versus f_PLL, or the bandwidth required, under the 1.9 kHz budget. In the papers read (Patra, Gong 2022 and 2023 and Gong's thesis, van Dijk 2019 in its published version, Ball et al. 2016 and the Lei 2024 review), none gives f_c,max versus f_PLL under the 1.9 kHz budget, checks the ≤ 1 MHz claim, or sweeps f_b. Lei cites van Dijk only as a pointer. For scale, Lei (p. 18) reports a 165–490 kHz corner for a SiGe VCO (Peng, ISSCC 2022; 13–17.5 GHz, bipolar, α not reported), against 1.3 MHz for Gong 2022 and 5.7 MHz for Patra. At f_PLL = 300 kHz a 165–490 kHz corner is tolerated for α ≤ 0.157–0.174 Hz²/Hz. Lei defines the corner as where flicker and thermal noise are equal and does not say it is the PN f⁻³ corner, so the match is likely, not certain. Ball et al. (main text and supplement) list flicker FM (p = −3) among the generic power-law regimes but give no corner tolerance. They apply a 1 Hz low-frequency cutoff that is negligible for their synthesizers, which need not hold for an oscillator with a MHz corner.

## 8. Files

Repository layout: src/budget.py (model, checkpoints, unit test), src/pll.py (linear PLL model), scripts/pll_comparison.py; scripts/digitize_patra_figs.py, scripts/fit_patra_data.py and scripts/fit_identifiability.py (the digitized data file is generated locally from the reader's own copy of Patra et al. and is not distributed); tests/test_budget.py; figures/budget_figures.png.
