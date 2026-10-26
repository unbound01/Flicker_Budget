# Phase 3: a generic linear PLL model and a three-case comparison

Code: `src/pll.py`, `scripts/pll_comparison.py`, `tests/test_pll.py`. Figures: `figures/pll_three_cases.png` (the comparison), `figures/pll_sensitivity.png`.

**Status:** scenario analysis with *generic* loop families. It is not a model of any specific PLL circuit, and its results are conditional on the assumptions below.

## Model

Output phase noise, single-sideband L (same convention as `budget.py`; FN² = 2∫f² L_out df, 8.3 kHz to 10 MHz in every case):

L_out(f) = |H_lp(f)|² · L_in(f) + |1 − H_lp(f)|² · L_osc(f)

- **L_osc**: the free-running 4 K oscillator at 6.3 GHz, either the two-term fit α/f²·(1 + f_c/f) (α = 0.137 Hz²/Hz, f_c = 4.12 MHz, free-running FN 3.28 kHz) or the digitized curve used directly above 150 kHz (free-running FN 3.27 kHz). No carrier scaling is applied: the 6.3 GHz spectrum is treated as the target oscillator.
- **L_in**: all input-side noise (reference × N², detector, divider, buffer) lumped into one flat floor, already referred to the output carrier, sent through one transfer function. Independent sources are assumed (no cross-spectra).
- **H_lp**: closed-loop transfer, with "bandwidth" meaning the closed-loop −3 dB frequency. Generic families: first-order; type-II with a PI zero (ζ = 1); second-order low-pass (ζ = 1.6). A steep 4th-order Butterworth loop is an extreme test. **Gong-matched case:** Gong's DAPLL is type-II (a V/I stage driving (R_Z + 1/sC_Z) ∥ 1/sC_P, thesis Fig. 6.6(a)) with ζ ≈ 1.6 and a 4 MHz bandwidth defined as the type-II closed-loop −3 dB point, so the closest generic family is type-II with ζ = 1.6, not the second-order low-pass (an earlier version of this note paired ζ = 1.6 with the low-pass family and called it Gong's). At ζ = 1.6 the two families put the natural frequency about 10× apart for the same bandwidth (BW = 3.51 f_n for type-II, 0.35 f_n for the low-pass). The C_P pole and the zero-order hold are not included.

### Floor scenarios (kept separate)

| Name | L_in | Like-for-like with |
|---|---|---|
| white-FM floor | α/f₃² (Patra's formula, derived assuming the corner is below the bandwidth) | Case 2 "fixed floor" |
| flicker-aware floor | L_osc(f₃) = α/f₃²·(1 + f_c/f₃), the full spectrum at the bandwidth | Case 2 "continuous floor" |
| absolute scenarios | flat −115, −122, −135, −140 dBc/Hz at 6.3 GHz | none; −115 is Patra's quoted number |
| perfect input | 0 | thought experiment only |

The **−122 dBc/Hz scenario** scales Gong 2023's −130 dBc/Hz at the 2.5 GHz divided output by 20·log₁₀(6.3/2.5). It assumes the same time-noise at both carriers and treats a measured closed-loop output floor as the input-referred floor, so it is an illustrative scenario, not a demonstrated floor of the modelled PLL. The same applies to the −141 dBc/Hz value obtained by scaling Gong's detector-only estimate (−137 dBc/Hz at 10 GHz).

The matched floors move with the bandwidth, so they are a comparative rule, not a fixed-noise design; they are unphysically high at low bandwidth.

Validation: with a brick-wall H_lp and the flicker-aware matched floor the code reproduces the closed-form budget (Case 2, continuous) to better than 0.5% (a unit test, not a result).

## Three cases, identical limits (fitted oscillator spectrum)

| Model | FN at 0.3 / 1 / 3 MHz (kHz) | Bandwidth for 1.9 kHz |
|---|---|---|
| Case 1: Patra, α/f² only (α = 0.154 from the measured PN at 10 MHz) | 1.73 / 1.69 / 1.57 | any (free-running already 1.75 kHz) |
| Case 2: brick wall, fixed floor α/f² | 2.58 / 2.27 / 1.89 | 2.91 MHz |
| Case 2: brick wall, continuous floor L_osc | 2.65 / 2.35 / 1.98 | 3.62 MHz |
| Case 3, white-FM floor: first-order / type-II / 2nd-order LP | 3.03 / 2.70 / 2.19 and 2.97 / 2.68 / 2.29 and 2.86 / 2.75 / 2.33 | 4.58 / 6.04 / 5.40 MHz |
| Case 3, flicker-aware floor: first-order / type-II / 2nd-order LP | 6.72 / 4.12 / 2.67 and 5.72 / 3.73 / 2.67 and 4.39 / 3.79 / 2.77 | 5.94 / 7.35 / 6.53 MHz |
| Case 3, Gong-matched (type-II, ζ = 1.6): white-FM / flicker-aware / perfect input | 3.00 / 2.69 / 2.23 and 6.26 / 3.93 / 2.66 and 2.61 / 2.28 / 1.85; at 4 MHz (Gong's bandwidth) 2.06 / 2.37 / 1.70 | 5.09 / 6.43 / 2.72 MHz |
| Case 3, scenario floor −122 dBc/Hz (any of the three loops) | 2.8 / 3.6–4.0 / 7.8–8.6 | never (minimum 2.75–2.80 kHz at 0.2–0.3 MHz) |
| Case 3, scenario floor −135 / −140 dBc/Hz | minima 2.33–2.48 / 2.05–2.19 kHz at 1.4 / 3.4–4.1 MHz | never |
| Case 3, floor −115 dBc/Hz (Patra's quoted number) | minimum 2.89–2.96 kHz at 0.1–0.2 MHz | never |
| Perfect input: Case 2 brick wall / first-order / type-II / 2nd-order LP | 2.57 / 2.25 / 1.81 for the brick wall | 2.52 / 2.27 / 3.51 / 3.31 MHz |

## Robustness checks

- **Direct digitized spectrum** (log-log interpolation above 150 kHz, fit below), bandwidth for 1.9 kHz [first-order / type-II / 2nd-order LP]: white-FM floor 4.40 / 5.75 / 5.16 MHz; measured-PN floor at the bandwidth 5.68 / 6.82 / 6.22 MHz; perfect input 2.20 / 3.40 / 3.14 MHz; scenario −122 dBc/Hz never. FN at 1 MHz: 2.7 kHz (white floor), 3.9–4.3 kHz (measured-PN floor), 3.6–4.0 kHz (−122), 2.2–2.5 kHz (perfect input).
- **Oscillator-fit windows** (16 windows; crossing computed per window, a sensitivity spread, not a confidence interval): white-FM floor 3.09–7.00 MHz (first-order), 4.11–8.19 MHz (type-II, 14 of 16 windows reach it), 3.90–7.84 MHz (2nd-order LP); flicker-aware floor 5.57–8.42 MHz (2nd-order LP); scenario −122 dBc/Hz never in any window.
- **Required floor** to meet 1.9 kHz: about −145 (type-II) and −143.6 (2nd-order LP) dBc/Hz at 5 MHz; no floor, even zero, works at 1 or 3 MHz for those two loops.
- **Matched-floor levels** at 0.3 / 1 / 3 MHz: white-FM −118.2 / −128.6 / −138.2 dBc/Hz; flicker-aware −106.5 / −121.5 / −134.4 dBc/Hz.

## Findings (conditional on these generic models)

1. **The brick-wall bandwidth (2.9-3.6 MHz) is an optimistic bound.** With Patra's own white-FM floor, generic loops need about 4.6-6.0 MHz (4.4-5.8 MHz from the direct spectrum), roughly 1.6-2.1× the brick-wall 2.91 MHz; with the flicker-aware floor 5.9-7.4 MHz, 1.6-2.0× the 3.62 MHz. Over the fit windows the spread is 3.1-8.4 MHz. For comparison, Gong's PLL uses 4 MHz with a 100 MHz reference. The Gong-matched loop (type-II, ζ = 1.6) needs 5.09 MHz with the white-FM floor and 6.43 MHz with the flicker-aware floor; at Gong's own 4 MHz it gives 2.06 kHz (white-FM floor), about 8% above the budget.
2. **At 1 MHz the target is not reached in any tested case:** 2.7 kHz with the white-FM floor, 3.6-4.1 kHz with the other floors, against 3.28 kHz free-running.
3. **Inside the −122 dBc/Hz scenario no bandwidth reaches 1.9 kHz**, and even −140 dBc/Hz gives 2.05-2.19 kHz at best. This is a statement about that scenario and these loops, not a hardware limit.
4. **A steeper roll-off does not help in the selected example.** The 4th-order Butterworth loop never gets below 2.85 kHz (perfect input) or 2.95 kHz (white-FM floor): the selected H_lp forces peaking in the complementary oscillator transfer 1 − H_lp. This does not show that every steep-roll-off architecture must fail.
5. **Case 1 cannot represent this oscillator:** α/f² only predicts the free-running oscillator already meets 1.9 kHz, contradicting the measured 3.4 kHz.
6. The perfect-input rows are thought experiments (all input noise set to zero, one chosen transfer function), not general lower bounds.

## What this does and does not show

- It shows that Patra's brick-wall reasoning does not carry over to generic smooth loops for this oscillator, and that under these models the achievable FN is limited on one side by in-band floor × bandwidth³ and on the other by the oscillator's flicker.
- It does **not** show that a PLL cannot meet 1.9 kHz. Patra measured no PLL; a specific loop with extra filtering, a shaped floor, separate injection points or injection locking could behave differently.

## Limitations and next step

Generic families, one aggregated floor through one transfer function, no cross-spectra, no circuit-specific loop-filter values, VCO gain, reference carrier or division ratio; the floor scaling to 6.3 GHz is my own conversion; the oscillator below 150 kHz is a fit extrapolation (matters little for bandwidths of 0.3 MHz and up); the fit carries the sensitivity discussed in `identifiability.md`; the budget is quoted at 6 GHz and the spectrum is at 6.3 GHz. The natural next step is one concrete loop (the phase-domain model of Gong's PLL, `docs/phase4_parameter_sheet.md`) with each noise source on its own transfer path, calibrated first against Gong's own measured 4.2 K profile.

## Corrections made after an external audit

- The "matched floor" earlier described as Patra's own floor was the flicker-aware floor L_osc(f₃); Patra's formula is the white-FM α/f₃². Both are now separate scenarios, and the earlier "5.9-7.4 MHz with Patra's floor assumption" is replaced by 4.6-6.0 MHz (white-FM) and 5.9-7.4 MHz (flicker-aware).
- "Realistic PLL" became "generic linear PLL model", and the −122 dBc/Hz value is labelled a scenario rather than a demonstrated floor.
- The oscillator-fit spread now computes the crossing per window instead of reading it from an envelope of curves.
- (16 Oct) Gong's loop is type-II, not a second-order low-pass: the ζ = 1.6 low-pass rows are kept as a generic family, and a Gong-matched type-II ζ = 1.6 row is added. The direct-spectrum and fit-window numbers for that row come from `scripts/pll_comparison.py` when run with the local data file.
