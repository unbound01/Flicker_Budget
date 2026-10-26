# Phase 4: parameter sheet for the architecture-specific PLL model (Gong's DAPLL)

Target architecture: J. Gong, E. Charbon, F. Sebastiano, M. Babaie, "A Cryo-CMOS PLL for Quantum Computing Applications", IEEE JSSC 58(5), 1362-1375, 2023
(analog charge-domain sub-sampling PLL, 40 nm, 10 GHz, 4.2 K). Sources: the **paper** (journal page numbers) and the author's **thesis** (Jiang Gong, TU Delft; chapter 6 "A Cryo-CMOS Dynamic-Amplifier-Based PLL", printed pp. 96-116 = PDF pp. 107-127; PDF page = printed page + 11).

Revised 16 Oct 2026 after reading the clean JSSC copy and the thesis: loop type corrected to type-II (§1), N = 108 at the Fig. 14 point, bandwidth definition, proxies marked as different chips, second consistency check (§3).

Status: **G** given in the source, **D** derived by me from given values (check before relying on it), **P** proxy from a different chip, **M** missing.

## 1. Architecture and noise paths

| Item | Content | Source | Status |
|---|---|---|---|
| Signal chain | input buffer (INBUF) → dynamic-amplifier phase detector (two-stage DA + resampling stage) → fully differential folded-cascode V/I stage → passive loop filter → LC-VCO → on-chip ÷4 (test divider, outside the loop) | paper p. 1368-1369, Fig. 8(a) | G |
| Frequency control | background frequency-tracking loop (FTL); manual coarse tuning within ±0.5·F_REF at power-up | paper p. 1369 | G |
| Linear phase-domain model | reference phase ×N → Σ → K_PD(s) → K_V/I → F_LPF(s) → 2π·K_VCO/s → ÷4; noise injected at: input/reference (INBUF), DAPD, V/I + loop filter, VCO | paper p. 1369 Fig. 8(b); thesis PDF p. 117 Fig. 6.6(b) | G |
| K_PD(s) | K_PD · (1 − e^(−s·T_REF))/(s·T_REF), a zero-order hold; the sinc response of the windowed-current integration is neglected | paper p. 1369, Eq. (25) | G |
| Loop type | **type-II with an extra pole.** The V/I stage drives a differential passive filter: a series R_Z–C_Z branch (drawn as 0.5R_Z + C_Z + 0.5R_Z) in parallel with C_P (2C_P per side), i.e. Z_LF = (R_Z + 1/sC_Z) ∥ 1/sC_P. With the VCO integrator this is a type-II loop with a zero and a third pole. The "type-I" on p. 1363 applies only to the simplified specification derivation ("for the sake of simplicity"). Same topology as the sister CSPLL, which the thesis calls type-II (Eq. 5.23, ζ = 0.5·R·C·ω_n) | thesis PDF p. 117, Fig. 6.6(a); thesis Ch. 5 Eq. (5.23); paper p. 1363 | G (schematic); ζ definition carried over from Ch. 5 is D |
| Division ratio N | integer-N with F_REF = 100 MHz (Table II), so N = F_VCO / 100 MHz: **N = 100 at 10 GHz (Fig. 13), N = 108 at 10.8 GHz (Fig. 14).** No divider in the loop; the ÷4 is a test divider | paper Table II p. 1372; p. 1371 | D |

## 2. Loop values that are given

| Parameter | Value | Source | Status |
|---|---|---|---|
| Reference frequency F_REF | 100 MHz (external sinusoid reshaped by the INBUF) | paper p. 1363, 1368 | G |
| PLL bandwidth bw_PLL | 4 MHz selected; larger would hurt stability and spurs | paper p. 1363 | G |
| Bandwidth definition | BW = √(2ζ² + 1 + √((2ζ² + 1)² + 1)) · ω_n/2π, the closed-loop −3 dB point of the type-II H(s); at ζ = 1.6, BW = 3.51 f_n, so 4 MHz gives f_n ≈ 1.14 MHz. Stated for the CSPLL (Ch. 5); assumed to carry over to the DAPLL | thesis Ch. 5 footnote (printed p. 77) | G for CSPLL, D for DAPLL |
| Damping factor | about 1.6 (lower damping suppresses VCO flicker better at the cost of stability) | paper p. 1369, 1371 | G |
| Phase-margin cost of C_RS | degradation within 14° for a 4 MHz bandwidth, with C_RS = C_S2/4 | paper p. 1369 | G |
| First-stage K_PD | about 0.6 V/rad; second stage deliberately low-gain (large C_S2); overall K_PD set digitally by C_S1, C_S2 | paper p. 1369 | G (overall value M) |
| K_PD model | K_PD = 4·G_M·A_VCO·sin(0.5·ω·T_P)/(ω·C_S); A_VCO = 0.45 V, T_P = 30 ps, (W/L)₁,₂ = 1.2 µm/40 nm; K_PD drops about 20% at 4.2 K | paper p. 1366-1367, Eq. (12), Fig. 5 | G |
| Carrier / tuning range | 10 GHz; 8.9-11.1 GHz at 300 K, 9.4-11.6 GHz at 4.2 K | paper p. 1370 | G |
| Temperature | 300 K and 4.2 K; loop parameters re-optimised at 4.2 K | paper p. 1370-1371 | G |

## 3. Noise sources that are given

| Source | Level | Reference plane | Source | Status |
|---|---|---|---|---|
| DAPD reset noise | 4KT/(C_S·F_REF) (low-frequency, before K_PD) | sampler output | paper p. 1367, Eq. (18) | G |
| DAPD (M₁,₂ thermal + shot + flicker) | in-band floor −128 dBc/Hz at 300 K, estimated −137 dBc/Hz at 4.2 K (shot-limited); flicker included in Eq. (22) but no device model | referred to a 10 GHz carrier | paper p. 1368 (Eq. 22), 1369 (levels) | G |
| INBUF | simulated floor −161 dBc/Hz at RT, limited to −167 dBc/Hz at 4.2 K (shot-limited); flicker dominates the in-band PN below 100 kHz | its own plane, not stated; I assume the 100 MHz reference | paper p. 1369, 1371 | G (plane D) |
| Measured in-band floor | −130 dBc/Hz at the 2.5 GHz divided output (300 K), limited by the INBUF | divided output (÷4) | paper p. 1370 | G |
| VCO | flicker corner "a few MHz" at 4.2 K, "partially mitigated by using a large DAPLL loop bandwidth"; free-running PN shown in Fig. 14(b) only as a curve | 10.8 GHz carrier (÷4 output at 2.7 GHz: add 12.04 dB to refer to 10.8 GHz) | paper p. 1370-1371, Fig. 14(b) | G (numbers only as a figure) |
| Measured output | PN profiles at 300 K / 4.2 K, rms jitter 75.3 fs / 36.7 fs (10 kHz-30 MHz, Table II), integrated PN −49.5 / −55.8 dBc | 10 GHz | paper p. 1370-1372, Figs. 13, 15, Table II | G |
| Measured noise attribution at 4.2 K | below 100 kHz: INBUF flicker; 100 kHz–4 MHz: filtered VCO flicker and INBUF shot noise; a lower ζ would suppress the VCO flicker better at the cost of stability | Fig. 14(b), 10.8 GHz | paper p. 1371 | G |

**Consistency check (mine):** referring the INBUF floor to the output with 20·log₁₀(N) = +40 dB gives −121 dBc/Hz (RT) at 10 GHz, i.e. −133 dBc/Hz at the 2.5 GHz output, against the measured −130 dBc/Hz. That is within 3 dB, which supports reading the INBUF figures as being at the 100 MHz reference and multiplied by N² in the loop. It is an inference, not a statement in the paper.

**Consistency check 2 (mine):** the integrated PN and jitter agree under the single-sideband convention σ_φ = √(2·10^(IPN/10)): −49.5 dBc gives 75.4 fs and −55.8 dBc gives 36.5 fs at 10 GHz, against 75.3 and 36.7 fs in Table II. So Gong's integrated PN uses the same single-sideband L as this repo (FN² = 2∫f²L df).

## 4. Missing, and how to fill it

| Missing | Why it matters | Proxy available | How to fill it |
|---|---|---|---|
| K_VCO | sets the loop gain and the VCO-noise high-pass | 50 MHz/V in an analysis example of the sister charge-sampling PLL (thesis printed p. 79, PDF p. 90). **Different chip.** The thesis Ch. 7 DAPLL (dual-core VCO, ÷8, 120 MHz/V) is also a different chip; do not use either as this chip's value | not needed if the loop is parametrised by (bw, ζ); otherwise ask the authors |
| R_Z, C_Z, C_P values, V/I transconductance (topology now known, see §1) | set the zero, the third pole and the roll-off of the reference path | sister chip (CSPLL, **different chip**): total loop-filter capacitance 20.7 pF, 4-bit tunable compensation resistor (thesis printed p. 83, PDF p. 94); C₁ ≈ 0.015·C (thesis Ch. 5) as a starting ratio for C_P/C_Z | parametrise by (bw, ζ, C_P/C_Z) plus the zero-order hold; fit the ratio to Fig. 14(b); or ask the authors |
| C_S1, C_S2 and the overall K_PD | closed-loop gain at the optimum setting | C_S = 30 fF used in the K_PD simulation (paper Fig. 5); C_RS = C_S2/4 | back out from bw = 4 MHz and ζ = 1.6 |
| Measured loop bandwidth and phase margin at the optimum C_S2 | the 4 MHz is the design target; Fig. 14(a) shows three C_S2 settings | none | read from Fig. 14(a) (digitize locally) |
| Numeric free-running VCO and per-source PN curves | calibrates the separate noise paths | none | digitize Fig. 14(b) (estimated contributions: reference, VCO, DAPD, total) and Fig. 13(b) from your own copy; keep the points out of the public repo |
| Reference-plane definition of the INBUF and DAPD floors | decides whether N² applies | consistency check above | ask the authors, or read thesis section 6.3 |
| Cryogenic flicker coefficients | the paper notes there is no flicker device model | none | take the measured curves instead of a model |

## 5. What this gives the Phase 4 model

- The **loop shape** can be fixed from given values: type-II with Z_LF = (R_Z + 1/sC_Z) ∥ 1/sC_P, bandwidth 4 MHz (type-II −3 dB definition), damping 1.6, zero-order hold at 100 MHz, N = 108 at the Fig. 14 point. The closest Phase 3 family is "type-II (PI zero), ζ = 1.6" (`pll.GONG_MATCHED`), not the second-order low-pass used there earlier. The C_P/C_Z ratio is the one free shape parameter.
- The **noise paths** are separable as in Fig. 8(b): reference/INBUF (×N²), DAPD, V/I plus loop filter, VCO, each through its own transfer function. The levels in section 3 are usable; the plane of the INBUF floor needs confirming.
- The **calibration target** is Gong's own measured and estimated profiles (Fig. 13(b), Fig. 14(b)), so the model can be checked against a real 4.2 K PLL before it is paired with Patra's oscillator.

## 6. A correction to Phase 3

The "−122 dBc/Hz Gong-scaled" scenario is the **300 K** measured floor (−130 dBc/Hz at 2.5 GHz) scaled to 6.3 GHz. For 4.2 K the INBUF estimate is 6 dB lower (−167 vs −161 dBc/Hz), which scales to about −131 dBc/Hz at 6.3 GHz (N = 63), and the DAPD-only estimate (−137 dBc/Hz at 10 GHz) scales to about −141 dBc/Hz. These fall between the −122 and −135 dBc/Hz curves already computed (minimum FN 2.8 and 2.4-2.5 kHz); I have not re-run the model at −131 dBc/Hz.
