# Is 1.9 kHz the right number for any operation? Three criteria on the same spectra

Code: `scripts/filter_criteria.py` (needs the locally generated `data/gong_fig13*.json` and `data/patra_digitized.json`), tests in `tests/test_filters.py`. Builds on `docs/a0_gong_measured.md`.

## Criteria

With S_f(f) = 2 f² L(f) the one-sided frequency-noise PSD (L single-sideband):

| Criterion | Definition | Pass level |
|---|---|---|
| **brick** (Patra) | σ² = ∫ S_f df over 8.3 kHz–10 MHz | σ ≤ 1.9 kHz |
| **gate** (Gong Eqs. 1–2, π rotation) | 1−F = ∫ S_f · \|H_LO\|² / f_R² df | 1−F ≤ 10⁻⁵ |
| **idle** (free evolution / Ramsey of length τ) | var_φ(τ) = (2πτ)² ∫ S_f · sinc²(πfτ) df | var_φ ≪ the qubit's own (2π · 1.9 kHz · τ)² |

The idle filter sees all frequency noise below about 1/(2τ) as quasi-static, so for a free-running oscillator var_φ grows as τ². Far above 1/τ it sees only the phase error at the start and end of the idle, so for a phase-locked LO var_φ saturates at 2σ_φ², with σ_φ² = 2∫L df the integrated phase noise. The 1.9 kHz figure is the qubit's own quasi-static frequency noise (T₂* = 1/(√2 π σ) = 118 µs; van Dijk Table II), so "σ_eff(τ) = √var_φ / (2πτ) versus 1.9 kHz" is the like-for-like use of Patra's number.

Self-checks (unit tests, not results): quasi-static limit; T₂* from 1.9 kHz; the white-PM limit 2σ_φ²; the gate filter's DC gain; the white-FM π-pulse value 9.08×10⁻⁷ already in `report.md`.

## Spectra

| Label | Spectrum | Carrier |
|---|---|---|
| G4 | Gong Fig. 13(b), measured PLL output, 4.2 K | 10 GHz |
| G300 | Gong Fig. 13(a), measured, 300 K | 10 GHz |
| P-FR | Patra 4 K oscillator, two-term fit, free running | 6.3 GHz |
| P-L1, P-L4 | the same oscillator in the Gong-matched loop (type-II, ζ = 1.6), BW 1 or 4 MHz, white-FM floor | 6.3 GHz |
| X-8.6 | **third chip:** Xiao et al., IEEE RFIC 2025 (Xidian), 28 nm fractional-N CP-PLL, integer mode, 4 K (19.4 fs) | 8.6 GHz |
| U-5.7 … U-6.9 | **fourth chip:** Peng et al., IEEE RFIT 2024 (UESTC), 65 nm sub-sampling PLL, integer-N, 4 K (~95 fs) | 5.7, 6.0, 6.3, 6.9 GHz |
| T-8.0, T-10.6, T-15.4 | **second chip:** Gong thesis Fig. 7.21, octave DAPLL with a dual-core triple-mode VCO (not the JSSC chip), measured at 4.2 K after ÷8 | 8.0, 10.6 and 15.4 GHz |

**Second-chip validation** (`scripts/digitize_gong_thesis_fig721.py`): every panel reproduces its own analyzer box over 10 kHz–30 MHz. Integrated noise is within 0.01 dB, PM within 0.2%, FM within 0.8%, and jitter matches (47.3, 40.6, 47.0 fs). Spot markers match within 0.17 dB, except one 10 kHz edge point at 0.54 dB.

Below Gong's first data point (10.5 kHz), three extrapolations down to 1 Hz: A, the measured slope (−8.8 dB/dec, near flicker PM, which Gong attributes to the INBUF); B, white FM; C, flicker FM below 1 kHz. **The results below are the same to two digits for A, B and C**: the idle filter cancels slow phase wander that is common to the start and end of the idle.

## Results

| Spectrum | brick FN | gate 1−F, f_R = 1 / 10 MHz | idle var_φ, τ = 1 µs / 120 µs [rad²] | idle relative to the qubit's own, 1 µs / 120 µs |
|---|---|---|---|---|
| G4 (measured, 4.2 K) | 8.86 kHz (**4.7× over**) | 5.5×10⁻⁶ / 1.3×10⁻⁶ (pass) | 1.4×10⁻⁵ / 1.4×10⁻⁵ | 0.1 / 7×10⁻⁶ |
| G300 (measured, 300 K) | 16.6 kHz (8.7× over) | 2.0×10⁻⁵ / 5.1×10⁻⁶ (fails at 1 MHz) | 4.6×10⁻⁵ / 4.6×10⁻⁵ | 0.3 / 2×10⁻⁵ |
| T-8.0 (2nd chip, 4.2 K) | 5.71 kHz (3.0× over) | 5.1×10⁻⁶ / 6.5×10⁻⁷ (pass) | 1.1×10⁻⁵ / 1.3×10⁻⁵ | 0.08 / 6×10⁻⁶ |
| T-10.6 (2nd chip, 4.2 K) | 5.56 kHz (2.9× over) | 7.3×10⁻⁶ / 5.8×10⁻⁷ (pass) | 1.4×10⁻⁵ / 1.8×10⁻⁵ | 0.1 / 9×10⁻⁶ |
| T-15.4 (2nd chip, 4.2 K) | 8.99 kHz (4.7× over) | **2.0×10⁻⁵ (fails)** / 1.6×10⁻⁶ | 4.0×10⁻⁵ / 4.8×10⁻⁵ | 0.3 / 2×10⁻⁵ |
| X-8.6 (3rd chip, 4 K) | 2.43 kHz (1.3× over) | 7.9×10⁻⁷ / 1.4×10⁻⁷ (pass, 13× margin) | — / 2.5×10⁻⁶ | — / 1×10⁻⁶ |
| U-5.7 (4th chip, 4 K) | 5.49 kHz (2.9×) | 9.9×10⁻⁶ (borderline) / 6.3×10⁻⁷ | — / 2.5×10⁻⁵ | — / 1×10⁻⁵ |
| U-6.0 | 6.57 kHz (3.5×) | **1.1×10⁻⁵ (fails)** / 8.7×10⁻⁷ | — / 2.7×10⁻⁵ | — / 1×10⁻⁵ |
| U-6.3 | 7.57 kHz (4.0×) | **1.3×10⁻⁵ (fails)** / 1.1×10⁻⁶ | — / 3.1×10⁻⁵ | — / 2×10⁻⁵ |
| U-6.9 | 9.17 kHz (4.8×) | **1.5×10⁻⁵ (fails)** / 1.6×10⁻⁶ | — / 3.3×10⁻⁵ | — / 2×10⁻⁵ |
| P-FR (free running) | 3.28 kHz (1.7× over) | 1.7×10⁻⁵ (fails) / 2.6×10⁻⁷ | 5.8×10⁻⁴ / **5.2** | **4 / 3** |
| P-L1 (1 MHz loop) | 2.69 kHz (1.4× over) | 2.0×10⁻⁶ / 1.5×10⁻⁷ (pass) | 8.4×10⁻⁶ / 9.1×10⁻⁶ | 0.06 / 4×10⁻⁶ |
| P-L4 (4 MHz loop) | 2.06 kHz (1.1× over) | 4.0×10⁻⁷ / 1.2×10⁻⁷ (pass) | 9.3×10⁻⁷ / 9.0×10⁻⁷ | 0.007 / 4×10⁻⁷ |

For the free-running oscillator the idle criterion gives σ_eff ≈ 3.0 kHz at τ = 120 µs, close to its brick value of 3.28 kHz. For every locked spectrum, σ_eff falls as 1/τ: 586 Hz at 1 µs and 4.9 Hz at 120 µs for G4.

**Robustness of the in-loop case to the floor** (Patra's oscillator, Gong-matched loop): with the white-FM, flicker-aware, −122 and −131 dBc/Hz floors, the gate infidelity at f_R = 1 MHz is (1.8–3.3)×10⁻⁶ at BW = 1 MHz and (0.4–6.8)×10⁻⁶ at BW = 4 MHz. Idle var_φ at 120 µs is (0.9–16)×10⁻⁶ rad². The brick integral over the same cases is 2.1–10.1 kHz. The tightest case is the −122 dBc/Hz floor at 4 MHz: gate 6.8×10⁻⁶ against the 10⁻⁵ target.

**Where each criterion looks** (G4, share of the criterion by offset band):

| Criterion | <10 kHz | 10–100 kHz | 0.1–1 MHz | 1–4 MHz | 4–10 MHz | >10 MHz |
|---|---|---|---|---|---|---|
| brick 8.3k–10M | 0% | 0% | 1% | 23% | **77%** | 0% |
| gate, f_R = 1 MHz | 0% | 0% | 12% | **66%** | 20% | 2% |
| idle, τ = 120 µs | 2% | 3% | 24% | **41%** | 27% | 4% |

**Second chip, gate infidelity versus rotation speed** (f_R = 0.5 / 1 / 2 / 4 / 10 MHz): 8.0 GHz 4.1, 5.1, 5.1, 2.3, 0.65 ×10⁻⁶; 10.6 GHz 6.0, 7.3, 5.1, 2.1, 0.58 ×10⁻⁶; 15.4 GHz 18, 20, 13, 5.3, 1.6 ×10⁻⁶. At 15.4 GHz the 99.999% target is met only for f_R ≥ 4 MHz.

**Third and fourth chips (independent groups).** Validation (`scripts/digitize_other_cryo_plls.py`): Xiao's 4 K trace reproduces its own analyzer box over 12 kHz–20 MHz (−62.57 vs −62.60 dBc, PM +0.2%, FM −0.2%, 19.46 vs 19.485 fs). Peng's four 4 K curves reproduce their printed PN at 100 kHz within 0.5 dB and jitter (1 kHz–100 MHz) within 1.8–3.7%; these figures are ~350 px bitmaps, so this chip is the least precise. Gate infidelity versus f_R (0.5 / 1 / 2 / 4 / 10 MHz): Xiao 0.83, 0.79, 0.61, 0.36, 0.14 ×10⁻⁶; Peng 5.7 GHz 11, 9.9, 4.9, 2.0, 0.63 ×10⁻⁶; 6.9 GHz 13, 15, 10, 5.3, 1.6 ×10⁻⁶.

**Across all nine 4 K points from four chips:** brick fails all nine (1.3–4.8×). The gate criterion at f_R = 1 MHz passes five (one borderline) and fails four; at f_R = 4 MHz all nine pass (≤ 5.3×10⁻⁶). Idle error is 2.5×10⁻⁶–4.8×10⁻⁵ rad² for all of them. Brick and gate infidelity are correlated (Spearman 0.72) but 8 of 36 pairs are ranked in opposite order.

**The brick metric does not rank PLLs correctly.** The JSSC chip at 10 GHz has a worse brick value than the thesis chip at 10.6 GHz (8.86 against 5.56 kHz) but a better gate infidelity at f_R = 1 MHz (5.5 against 7.3 ×10⁻⁶).

## What this shows (conditional on the limits below)

1. **1.9 kHz is the right number for one thing: quasi-static frequency noise during free evolution.** It is a valid test of a free-running oscillator. There the idle criterion and the brick integral roughly agree (3.0 against 3.3 kHz), and both say Patra's free-running oscillator would roughly triple the qubit's dephasing.
2. **For a phase-locked LO it is the wrong metric.** The loop turns low-frequency frequency noise into bounded phase jitter. Idle error then becomes a τ-independent phase-jitter term, about 2σ_φ² ≈ 1.4×10⁻⁵ rad² for Gong at 4.2 K, which is set by the integrated jitter and not by an FN integral. The gate error is dominated by 1–4 MHz, while the brick integral is dominated by 4–10 MHz.
3. **Under both physical criteria, Patra's ≤1 MHz claim is supported for his own 4 K oscillator** in a Gong-type loop, for every floor scenario tested. Only his brick metric fails it, by 1.3–2×. This reverses the reading of the earlier result "≤1 MHz is not reproduced". That result stands, but only for the brick metric. The f_c tolerances in `report.md` answer "how much flicker does the brick budget tolerate", which is not the question that matters for a locked LO.
4. **Four measured cryogenic PLLs from three groups (nine 4 K operating points) fail the brick budget at every point, by 1.3–4.8×, including the lowest-jitter cryo PLL in our set (19.4 fs), which meets the gate criterion with a 13× margin.** The details for the first two chips follow. Under the gate criterion at f_R = 1 MHz, three of the four points pass and one (15.4 GHz) fails. The idle error is about 10⁻⁵ rad² at all four, under 3×10⁻³ of the qubit's own dephasing for τ ≥ 10 µs. The brick budget therefore does not discriminate between these PLLs, while the gate criterion does. At 300 K the JSSC PLL fails the gate criterion at f_R = 1 MHz.
5. For a locked LO, the useful specifications are **integrated jitter** (idle) and the **gate-filtered FN** (driven), not a broadband FN integral.

## Limits

- **Dynamical decoupling is not covered.** A CPMG or XY sequence with pulse spacing t_p has a band-pass filter peaked near 1/(2t_p). For t_p = 0.1–1 µs that falls at 0.5–5 MHz, which is exactly the region the brick integral weights. This is the one operation for which MHz frequency noise could matter, and it is the next check.
- Single gate, worst-case π-rotation filter. No two-qubit gates, no accumulation over long sequences (coherent errors may add faster than these incoherent estimates), and a Gaussian, stationary noise model.
- The Gong spectra below 10.5 kHz are extrapolated (the three scenarios agree). Spectra above 30 MHz are extrapolated with the measured end slope.
- **Gate integration band.** Gate infidelity is integrated over each measured band (to 30 MHz for Gong, about 95 MHz for Xiao and Peng), which reproduces Gong's stated fidelity. The filter still weights far-out phase noise roughly one-to-one: a flat floor at the last measured level out to 1 GHz raises 1 − F by 4–36% and moves only Peng 5.7 GHz across 1e-5 (9.9 → 10.4 ×10⁻⁶); out to 10 GHz all Gong points would fail while Xiao's still passes (2.8×10⁻⁶). The verdicts assume the drive chain band-limits the LO noise (Ball 2016); alternative-integrator recomputation: `scripts/alternative_integrator.py`.
- The Patra in-loop spectra use a generic type-II loop and assumed floors (see `docs/pll_model.md`), not a validated circuit model.
- The carriers differ: Gong at 10 GHz, Patra at 6.3 GHz. The comparison within each spectrum is unaffected.
- The 4.2 K Gong trace does not reproduce its own integrated box (`a0_gong_measured.md`). Scaling it down by that mismatch makes the margins larger, not smaller.

## Echo and worst-case idle (added for draft4)
- `echo_var(L, tau)`: Hahn-echo phase variance, weight 16 sin^4(pi f tau/2) on S_phi. Far-offset weight averages 6 (Ramsey: 2), so locked LOs get ~3x worse; free-running oscillators get ~10x better.
- `idle_infidelity(var, worst_case)`: average (1-e^{-v/2})/3 or worst case (equatorial state) (1-e^{-v/2})/2.
- At 120 us: Patra free-running 1-F 0.31 (echo 0.067); Gong 2022 at 6.3 GHz 0.12 (echo 0.013), at 4.65 GHz 0.073 (echo 0.007). A 1e-5 average idle target needs sigma_eff ~ 10 Hz.
