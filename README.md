<!-- DRAFT written by AI as a starting point. Rewrite in my own voice before anyone reads it (see SPEC.md Section 9). Read it out loud; change anything that doesn't sound like me. -->

# Predicting the LCST of NIPAM copolymers

PNIPAM dissolves in cold water and turns cloudy at about 32 C. Copolymerising NIPAM with a second monomer moves that temperature, and the smart nanofibers made in Mitsuhiro Ebara's group at NIMS depend on this: they use NIPAM with N-hydroxymethylacrylamide (HMAAm) to push the switch up to the mild-hyperthermia range, around 42 to 45 C.

I wanted to see how far a simple model gets at predicting that temperature from the comonomer's structure and how much of it is in the chain. I collected published LCSTs for NIPAM copolymers, described each comonomer with a few descriptors, and checked the model on comonomers it had never seen, including Ebara's own HMAAm copolymers.

## Main result

Four fraction-weighted features (how much the comonomer changes logP, polar surface area and H-bond donors compared with NIPAM) in a ridge regression predicted the LCST of held-out comonomers with a mean error of 4.3 C. Predicting 32 C for everything gives 8.9 C, and a straight line in mole fraction gives 6.6 C.

On Ebara's HMAAm copolymers, with every HMAAm row removed from training, the error was 2.6 C across 7 measurements, and all 7 fell inside the 90% interval. The model runs about 3 to 4 C high near 20 mol% HMAAm.

The 4215 alvaDesc descriptors did not help. They looked fine under random cross-validation and fell apart on new comonomers.

![Held-out predictions](results/figures/parity_logo.png)

## How to run

Python 3.13 (the spec said 3.11, which isn't on this machine; noted in the SPEC changelog).

```bash
pip install -r requirements.txt
bash run.sh
```

`run.sh` rebuilds everything from `data/raw/` in about 7 minutes. Stage 4 is the slow one.

## What each stage does

| Script | Does | Writes |
|---|---|---|
| `src/s1_comonomers.py` | checks the dataset schema, parses every comonomer SMILES, draws them | `data/processed/comonomers.csv`, `results/figures/comonomer_structures.png` |
| `src/s2_descriptors.py` | RDKit descriptors, fraction-weighted features, reads the alvaDesc export | `data/processed/model_table.csv`, `descriptors_rdkit.csv`, `descriptors_alvadesc.csv` |
| `src/s3_clean.py` | drops missing, constant, near-constant and correlated descriptors (4226 to 470) | `data/processed/features_clean.csv`, `results/feature_cleaning.csv` |
| `src/s4_model.py` | baselines, ridge, random forest, boosting; random 5-fold, grouped 5-fold, leave-one-comonomer-out | `results/metrics.csv`, `results/predictions_logo.csv` |
| `src/s5_trends.py` | logP trend at fixed fraction, composition series, coefficients | `results/trend*.csv`, `results/coefficients_ridge_physics.csv`, figures |
| `src/s6_ebara_predict.py` | Ebara HMAAm check with HMAAm held out, new copolymer predictions, applicability domain | `results/ebara_check.csv`, `results/new_predictions.csv`, `results/applicability_domain.csv`, figures |

## Folders

```
data/raw/          lcst.csv (the dataset), alvadesc_export.txt (raw descriptor export)
data/processed/    comonomers, descriptors, model table, cleaned features
data/sources.md    every source, what was taken from it, and what was left out and why
data/pdfs/         paper PDFs used for extraction (git-ignored, not in the repo)
results/           metrics, predictions, trend tables, figures/
notes/             research_note.md
SPEC.md            the plan I followed, with the changelog
```

## About the data

- 113 training rows from 22 sources, covering 22 comonomers and 16 PNIPAM reference rows. 7 more rows from four Ebara group papers are kept in the same file with `role = ebara_check` and never used for training.
- Most rows come from tables. 9 were read off figures (`read_method = digitized_from_figure`), good to about 1 C.
- 19 rows are values another paper compiled from the literature (Hayakawa et al. 2023, ESI). The method and whether the fraction is feed or measured composition aren't given there, so these carry a `secondary` note.
- LCST is measured in different ways across papers: 50% transmittance, first visible turbidity, absorbance thresholds, DSC. The method is recorded for every row. It adds a few degrees of noise that the model can't remove.
- Charged comonomers (acrylic acid and two amine monomers) change LCST with pH, so the headline numbers leave them out. Results with them included are in `metrics.csv` under `subset = all`.

## Known limits

- 22 comonomers is not many, and 18 of them come from a single paper each.
- Crippen logP gives N-n-propylacrylamide the same value as NIPAM, although its homopolymer switches about 10 C lower. Simple 2D descriptors miss that kind of difference.
- The TPSA term in the model has a negative coefficient. It is strongly correlated with the logP term (r = -0.74), so the individual coefficients shouldn't be read as separate physical effects.
