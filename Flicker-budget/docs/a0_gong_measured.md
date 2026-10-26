# A0: Patra's 1.9 kHz budget versus a measured cryogenic PLL (Gong et al., JSSC 2023)

Code: `scripts/digitize_gong_fig13b.py`, `scripts/a0_gong_measured.py`, `scripts/digitize_gong_fig13a_check.py`. Data: generated locally from your own copy of the paper (`data/gong_fig13*.json`, git-ignored).

**Question.** Does a real, measured 4.2 K cryo-CMOS PLL meet Patra's 1.9 kHz rms frequency-noise (FN) budget over 8.3 kHz to 10 MHz, and how does that compare with the gate-fidelity criterion the same PLL was designed for? No loop model and no fitted spectrum are used here.

**Source.** Gong et al., Fig. 13(b): analyzer screenshot, PLL locked at ≈10 GHz, 4.2 K, measured after the on-chip ÷4 (2.500018 GHz). The screenshot prints spot-noise markers and the instrument's own integrated results over 10 kHz to 30 MHz. Fig. 13(a) (300 K) is used as a second validation case.

## Validation

| Check | Digitized | Reference | Result |
|---|---|---|---|
| 4.2 K spot noise at 100 kHz, 1 MHz, 10 MHz (2.5 GHz) | −132.72, −133.83, −149.01 dBc/Hz | −132.74, −133.87, −149.01 | within 0.04 dB |
| 4.2 K spot noise at 10 kHz (trace edge) | −127.25 | −127.03 | 0.22 dB |
| 300 K spot noise, 4 markers | | | within 0.25 dB |
| **300 K integrated, 10 kHz–30 MHz** | −61.55 dBc, PM 1183 µrad, FM 5982 Hz, 75.3 fs | −61.55 dBc, 1180 µrad, 5949 Hz, 75.36 fs | **within 0.5%** |
| 4.2 K integrated, 10 kHz–30 MHz | −66.83 dBc, PM 644 µrad, FM 2930 Hz, 41.0 fs | −67.80 dBc, 576 µrad, 2308 Hz, 36.66 fs | **+0.97 dB, FM +27%: not reproduced** |
| Gong's fidelity (Eqs. 1–2, θ = π), trace only | 99.99945% (f_R = 1 MHz), 99.99987% (10 MHz) | 99.9994%, 99.9998% | 1 MHz consistent; 10 MHz off |
| Same, plus the measured −69 dBc reference spur at 100 MHz | 99.99942%, 99.99985% | 99.9994%, 99.9998% | both consistent (10 MHz only if Gong truncated) |

Reading: the method and the conventions (SSB L; int noise = ∫L df; PM = √(2∫L df); FM = √(2∫f²L df)) are validated by the 300 K case. At 4.2 K the printed integrated values do not match the plotted trace, although the markers do; the cause cannot be determined from the figure. The FN results below are therefore given as a range: nominal (from the trace) and scaled down by the instrument's FM ratio (a crude lower bound). The spur reconciliation of Gong's fidelity is an inference: it assumes −69 dBc per sideband at 100 MHz and the worst-case |H_LO|².

## Results (4.2 K)

**Integrated FN, 8.3 kHz–10 MHz (target 1.9 kHz rms):**

| Reference plane | FN (trace) | FN (scaled to instrument) | vs 1.9 kHz |
|---|---|---|---|
| 2.5 GHz test output, as measured | 2.22 kHz | 1.74 kHz | 1.17× / 0.92× |
| 10 GHz PLL output (×4), the LO itself | **8.86 kHz** (8.56–9.17 for ±0.3 dB) | **6.98 kHz** | **4.7× / 3.7× over** |
| 6 GHz carrier with the same time noise (×6/2.5) | 5.32 kHz | 4.19 kHz | 2.8× / 2.2× over |

The 10 GHz plane is the relevant one for a 10 GHz spin qubit driven directly by this LO; dividing the carrier by N lowers FN by N, so the 2.5 GHz value does not describe the LO. At 300 K the same integral is 16.6 kHz.

**Where FN² comes from (10 GHz):** 8.3–100 kHz 0%, 0.1–1 MHz 1%, 1–4 MHz 23%, **4–10 MHz 77%**. Cutting the band: 0.75 kHz to 1 MHz, 1.75 kHz to 2 MHz, 4.4 kHz to 4.05 MHz. The budget is set almost entirely by the region just above Gong's 4 MHz loop bandwidth, where the loop hands over to the VCO (its flicker region at 4.2 K) and the f² weighting is largest.

**Gate fidelity of the same spectrum (Gong Eqs. 1–2, π rotation, 8.3 kHz–30 MHz, no spur):** 1−F = 5.0, 5.5, 6.1, 5.0 and 1.3 ×10⁻⁶ for f_R = 0.5, 1, 2, 4 and 10 MHz. The 99.999% target is met at every rotation speed. At f_R = 1 MHz only 20% of the infidelity comes from 4–10 MHz offsets, against 77% of FN².

## What this shows

1. **A measured 4.2 K cryo-CMOS PLL that meets a 99.999% single-gate fidelity target exceeds Patra's 1.9 kHz budget by about 4–5× at its 10 GHz output.** This does not depend on any fit, loop model or floor assumption; it depends on the digitized curve, with the 4.2 K box mismatch bounding the error (3.7–4.7×).
2. **The two criteria disagree because they weight offsets differently.** The brick-wall FN integral weights every offset by f² up to 10 MHz, so it is dominated by 4–10 MHz. The gate filter |H_LO|² is a low-pass of about 1.9 f_R (Gong p. 1363), so for f_R ≈ 1 MHz most of that region is rejected.
3. **This reframes the project's question, not just its numbers.** The 1.9 kHz figure is the qubit's own quasi-static frequency noise for T₂* = 120 µs (van Dijk Table II). It is a "do not degrade free evolution" criterion. Patra applies it to a 10 MHz broadband integral, where it acts as a much stricter driven-gate criterion than the gate filter itself. The corner tolerances in `report.md` are correct under that budget, but the budget is conservative for driven gates by a large factor and of unclear relevance for free evolution, whose filter (low-pass at ~1/T) weights offsets below the 10 kHz the data reach.

## Limits

- One measured curve, one operating point (10 GHz, optimum C_S2, 4.2 K). 300 K is a second point.
- 4.2 K integrated box not reproduced from its own trace (+0.97 dB); results are given as a range.
- 8.3–10 kHz is extrapolated with the local slope (−8.8 dB/dec); it carries under 0.1% of FN².
- Fidelity uses Gong's worst-case π-rotation filter for one gate; it is not a full error budget (no idle/Ramsey filter, no echo, no multi-gate accumulation).
- The 10 GHz plane assumes the LO is used at the VCO frequency; Patra's budget was quoted for 6 GHz.

## Next

- Replace the 1.9 kHz brick-wall integral by filter-function criteria in the analysis: Gong's gate filter (driven) and a Ramsey/idle filter (free evolution), evaluated on the same spectra. This is the flicker-aware fidelity budget Gong's footnote 2 leaves to simulation.
- Then Experiment A (loop reproduction from Fig. 14(b)) to sweep bandwidth and damping under those criteria.
