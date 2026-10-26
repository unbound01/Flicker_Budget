# Is a broadband frequency-noise budget the right specification for cryo-CMOS qubit LOs?

Analysis code for the paper by Krishna Gupta (Rajiv Gandhi Institute of Petroleum Technology, India):
*"Is a Broadband Frequency-Noise Budget the Right Specification for Cryo-CMOS Qubit Local Oscillators?
Reanalysis of Measured 4-K PLLs"* (preprint: [TechRxiv link to be added]).

The widely used budget of Patra et al. (JSSC 2018) asks a qubit LO for less than 1.9 kHz rms frequency noise
over 8.3 kHz-10 MHz. This repository

- digitizes the measured 4-K phase noise of four published cryo-CMOS PLLs (nine operating points) and validates
  each spectrum against the integrated noise, jitter or spot values printed by its authors;
- evaluates the same spectra with the brick-wall budget and with filter-function criteria for a single-qubit
  gate (Gong et al., JSSC 2023) and for an idle period (Ramsey and Hahn echo);
- extends the budget in closed form to oscillators with a 1/f³ flicker region (Appendix B of the paper).

## Reproduce the paper

```bash
pip install -r requirements.txt
# 1. regenerate the digitized data from your own copies of the source papers (not redistributed; see data/README.md)
python scripts/data_checksums.py      # checks the regenerated data against data/CHECKSUMS.sha256
# 2. tables, numbers and figures
python scripts/table3.py              # Table II rows + paper/generated/numbers.json (every number in Sections III-V)
python paper/make_figs.py             # all figures, into paper/figs/
python scripts/gate_two_level.py      # independent check of the gate criterion and the 2/3 state-average factor
python -m pytest                      # 27 tests (data-dependent tests are skipped if data/ is empty)
```

Appendix B numbers: `src/budget.py` (closed form, tolerable corner), `scripts/fit_patra_data.py` (fit),
`scripts/fit_identifiability.py` (fit windows, Monte Carlo), `scripts/pll_comparison.py` (loop models).

> **Status:** preprint, not yet peer reviewed. Spectra are digitized from published figures
> (about ±0.3 dB per point for vector figures, ±0.5 dB for bitmaps).

## Flicker-extension results (Appendix B; rows marked † need your own copy of the Patra paper; see `data/README.md`)

| Quantity | Value |
|---|---|
| α_max (f_PLL = 300 kHz, 8.3 kHz to 10 MHz, 1.9 kHz rms) | **0.1842 Hz²/Hz** (Patra: < 0.2) |
| In-band PN / PN at 10 MHz at α_max | −116.9 / −147.3 dBc/Hz (Patra: −115 / −147; the 2 dB in-band gap is unexplained but changes α_max by < 0.6%) |
| Corner tolerance f_c,max at α = 0.17 / 0.15 / 0.10 (fixed floor) | 233 / 637 / 2353 kHz (continuous floor: 213 / 582 / 2149 kHz) |
| Upper band limit f_b = 1 / 2 / 5 / 10 MHz | α_max = 2.26 / 1.00 / 0.38 / 0.18 Hz²/Hz (equivalent f_R = 0.41 / 0.81 / 2.03 / 4.05 MHz) |
| Gong 2022 oscillator (β = 0.041, f_C = 1.3 MHz) | FN = 1.09 kHz; f_c,max ≈ 9 MHz (about 7× margin) |
| † Fit to Patra's 4 K data (Fig. 22) | α ≈ 0.14 Hz²/Hz, f_c ≈ 4.1 MHz for the primary range, but **strongly range-dependent**: 0.07–0.24 Hz²/Hz and 2.0–10.5 MHz over 16 fit windows (a sensitivity spread, not a confidence interval), while K = α·f_c varies only 0.46–0.76×10⁶ Hz²; free-running FN 3.1–3.5 kHz vs Patra's 3.4 kHz |
| † Same oscillator through a PLL | brick-wall model: 2.11–2.56 kHz at 1 MHz in all 16 fit windows; the primary fit needs 2.9 MHz (fixed floor) to 3.6 MHz (continuous floor). **Generic linear-loop model (docs/pll_model.md): 4.6–6.0 MHz with Patra's white-FM floor α/f² (4.4–5.8 MHz from the direct digitized spectrum), 5.9–7.4 MHz with a flicker-aware floor; never reached under an illustrative scaled −122 dBc/Hz floor scenario** |

**Main findings**

1. Eq. 7 depends only on α and K = α·f_c, so f_c,max(α) is a hyperbola in (α, f_c).
2. For the digitized 4 K oscillator at 6.3 GHz, Patra's conclusion that a PLL of ≤ 1 MHz "can easily" reach 1.9 kHz is **not reproduced** under the extended spectrum and the tested PLL/floor models. This is a conditional model comparison, not a claim that Patra is wrong or that no PLL can meet the budget.
3. The data constrain K = α·f_c and the PLL-predicted FN far better than α and f_c separately (docs/identifiability.md). Exploratory only: K may rise on cooling for this chip, but the 165 K and 300 K fits are weakly constrained.
4. The brick-wall PLL is an optimistic bound: generic first-order, type-II and second-order loops need roughly 1.6-2× the bandwidth (4.4–7.4 MHz depending on floor assumption and loop family; 3.1–8.4 MHz over the fit windows). At 1 MHz these models give 2.7–4.1 kHz against 3.28 kHz free-running. The family closest to Gong's DAPLL (type-II, ζ = 1.6) needs 5.1–6.4 MHz and gives 2.06 kHz at Gong's 4 MHz (docs/pll_model.md).
5. The 1.9 kHz figure is the qubit's own frequency noise (van Dijk 2019, Table II), not a
   fidelity-derived oscillator limit.

Full write-up: [`docs/report.md`](docs/report.md). Baseline and data validation: [`docs/validation.md`](docs/validation.md). Parameter identifiability: [`docs/identifiability.md`](docs/identifiability.md). Generic PLL model: [`docs/pll_model.md`](docs/pll_model.md), with the three-case figure `figures/pll_three_cases.png`. Architecture-specific next step: [`docs/phase4_parameter_sheet.md`](docs/phase4_parameter_sheet.md).

## Quick start

```bash
pip install -r requirements.txt
python src/budget.py                 # baseline checkpoints + figures/budget_figures.png
python tests/test_budget.py          # 9 model tests (or: pytest); 4 data tests are skipped until step below
python tests/test_pll.py             # 6 PLL-model tests (no data needed)

# optional, needs YOUR copy of the Patra et al. PDF (not included):
python scripts/digitize_patra_figs.py --pdf PATH/TO/CryoCMOS_BPatra.pdf --out data/patra_digitized.json
python scripts/fit_patra_data.py     # fits + figures/patra_fit.png (git-ignored)
python tests/test_budget.py          # now 13 tests
```

### Reproduce the manuscript (needs the regenerated data in `data/`; check with `python scripts/data_checksums.py`)

```bash
python scripts/table3.py      # paper/generated/table3_rows.tex (Table III body) + paper/generated/numbers.json
python paper/make_figs.py     # paper/figs/fig1..fig7 (.pdf/.png) + fig_numbers.json
```
`numbers.json` holds every value quoted in Sections V-C to VII of the paper: both carriers, gate at 0.5-10 MHz,
Ramsey/echo and average/worst-case idle, the 1 GHz / 10 GHz floor sensitivities, the share of gate error above 3 f_R,
the free-running cases, the Spearman/discordant-pair counts and the jitter limits.

Tested with Python 3.12, numpy 2.4, scipy 1.17, matplotlib 3.10.

## Repository layout

```
src/pll.py                       linear phase-domain PLL model (loop transfer functions, integrated FN)
src/budget.py                    model: FN budget, closed form, f_c,max, alpha_max; checkpoint printout
scripts/pll_comparison.py         Phase 3: FN versus PLL bandwidth for the three models (figures/pll_*.png)
scripts/fit_identifiability.py   Phase 2: fit ranges, residuals, uncertainty, Monte Carlo (figure is git-ignored)
scripts/fit_patra_data.py        fit (alpha, f_c) to Patra Figs 22/24/25; 3.4 kHz and PLL checks
scripts/digitize_patra_figs.py   re-extract the data from YOUR copy of the paper (optional)
data/README.md                   how to generate patra_digitized.json locally (file not distributed)
tests/test_budget.py             regression tests (checkpoint values, units; fit tests need the local data)
docs/report.md                   the report
docs/validation.md               Phase 1 baseline and digitization checks
docs/identifiability.md          Phase 2 fit-range, uncertainty and trade-off study
docs/pll_model.md                Phase 3 generic PLL transfer-function model, three-case comparison
docs/phase4_parameter_sheet.md   Phase 4 parameter sheet for Gong's DAPLL (given / missing / proxy)
docs/a0_gong_measured.md         A0: Patra's 1.9 kHz budget vs Gong's measured 4.2 K PLL (Fig. 13), with gate fidelity
docs/filter_criteria.md          brick (1.9 kHz) vs gate-fidelity vs idle/Ramsey criteria on the same spectra
scripts/filter_criteria.py       the three-criteria comparison (tests/test_filters.py)
scripts/digitize_other_cryo_plls.py   third and fourth chips: Xiao RFIC 2025 (8.6 GHz) and Peng RFIT 2024 (5.7-6.9 GHz), 4 K, validated vs printed values
scripts/digitize_gong_thesis_fig721.py  second chip: Gong thesis Fig. 7.21 (4.2 K, 8.0/10.6/15.4 GHz), validated vs its analyzer boxes
scripts/digitize_gong_fig13b.py  A0 digitizer, Gong Fig. 13(b) (4.2 K); fig13a_check for 300 K; a0_gong_measured.py analysis
scripts/table3.py                Table III rows + numbers.json for the manuscript (tests/test_table3.py)
paper/make_figs.py               manuscript figures
figures/                         generated figures
```

## Model in one paragraph

Oscillator L(Δf) = α/Δf²·(1 + f_c/Δf) (single-sideband). A brick-wall PLL at f_PLL replaces the
oscillator noise below f_PLL with a flat floor, either α/f_PLL² (*fixed*) or α(1 + f_c/f_PLL)/f_PLL²
(*continuous*). FN² = 2∫Δf² L dΔf from 8.3 kHz to 10 MHz. Units are Hz throughout; α = L·Δf²,
and the one-sided FN spectrum is 2α. The conventions were cross-checked against van Dijk 2019
(Table II) and Ball et al. 2016 (SI Eqs. 80–83, 112).

## Data provenance and limitations

- The digitized points are **not included**. They are derived from figures of a paper copyrighted by IEEE,
  and IEEE says that reproducing its graphics and tables normally requires permission; I found no IEEE
  statement covering digitized data points, so the repo errs on the safe side (not legal advice).
  Generate the file locally from your own copy with `scripts/digitize_patra_figs.py`; it is git-ignored.
  The overlay images the script can write, and `figures/patra_fit.png`, redraw the paper's curves: do not publish them.
- `scripts/digitize_patra_figs.py` needs poppler's `pdftoppm` and OpenCV. Its pixel calibrations are specific
  to the 13-page IEEE Xplore PDF (pages 9 and 10 at 300 dpi).
- Results marked † in the table above and in the report are my own derived numbers (fit parameters, FN values)
  with the source cited; the underlying points are about ±0.3 dB.
- Main limitations: brick-wall PLL model; the data start at 100 kHz (about 31% of the 4 K
  free-running FN² is extrapolated, but the PLL results use only data above 300 kHz); the fixed
  floor is discontinuous at f_PLL (hence the continuous variant); 1.9 kHz is a heuristic budget, not
  a measurement or a fidelity limit; fits at 165 K and 300 K are weakly constrained. See report §6.
- The closed-form-versus-quad test checks the algebra; it is not a scientific result.

## References

1. B. Patra et al., "Cryo-CMOS circuits and systems for quantum computing applications," *IEEE J. Solid-State Circuits* 53(1), 309–321, 2018.
2. J. Gong et al., "A cryo-CMOS oscillator with an automatic common-mode resonance calibration for quantum computing applications," *IEEE Trans. Circuits Syst. I* 69(12), 4810–4822, 2022.
3. J. Gong et al., "A cryo-CMOS PLL for quantum computing applications," *IEEE J. Solid-State Circuits* 58(5), 1362–1375, 2023.
4. J. P. G. van Dijk et al., "Impact of classical control electronics on qubit fidelity," *Phys. Rev. Applied* 12, 044054, 2019 (arXiv:1803.06176).
5. H. Ball, W. D. Oliver, M. J. Biercuk, "The role of master clock stability in quantum information processing," *npj Quantum Inf.* 2, 16033, 2016 (arXiv:1602.04551).
6. L. Lei et al., "Cryogenic control and readout integrated circuits for solid-state quantum computing," arXiv:2410.15895.

## License

Code: MIT, © 2026 Krishna Gupta (see `LICENSE`). Digitized data: not distributed, see above.
