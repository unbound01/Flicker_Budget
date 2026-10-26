# data/

`patra_digitized.json` is **not distributed** (it is git-ignored). It contains points digitized from
figures 22, 24 and 25 of Patra et al., IEEE JSSC 53(1), 2018, which is copyrighted by IEEE.

Generate it from your own legitimately obtained copy of the paper:

```bash
python scripts/digitize_patra_figs.py --pdf PATH/TO/CryoCMOS_BPatra.pdf --out data/patra_digitized.json
```

(needs poppler's `pdftoppm` and OpenCV; calibrations are for the 13-page IEEE Xplore PDF).
Then `scripts/fit_patra_data.py` and the four data-dependent tests will run.
Without the file, those tests are skipped and the model tests still run.
