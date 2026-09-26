# SPEC: Predicting the LCST of thermoresponsive copolymers for Ebara's smart nanofibers

Status: draft v1
Owner: (you)
Last updated: 2026-09-24

---

## 0. How to use this spec

Same rules as my other specs:

1. The spec is the source of truth. Build stage by stage.
2. Don't start a stage until the one before passes its "Done when" checks.
3. If reality differs, change the spec and log it in the Changelog.
4. Small stages. Split anything that feels big.

Read Section 9 (writing rules) before writing any README, comments, note or email.

---

## 1. What I'm building and why

Mitsuhiro Ebara's group at NIMS makes "smart" polymers that switch behaviour with temperature, and spins them into nanofiber meshes for drug release and hyperthermia treatments. Most of these are built on poly(N-isopropylacrylamide), PNIPAM. PNIPAM is dissolved in cold water and collapses (turns cloudy, pushes out water) above a transition temperature called the LCST, around 32 C. Copolymerising NIPAM with another monomer moves that temperature: a more hydrophilic comonomer usually raises it, a more hydrophobic one lowers it. Ebara's nanofibers use copolymers of this kind (for example NIPAM with a crosslinkable comonomer such as N-hydroxymethylacrylamide; check exact systems on his SAMURAI publication list), and tuning the switch temperature to just above body temperature is the whole design.

Goal: build a dataset of NIPAM copolymers with measured LCSTs from the literature, describe each by its comonomer (descriptors) and composition, and fit a model that predicts LCST. Check it reproduces the known hydrophilic-up / hydrophobic-down trend and the composition dependence, check it on Ebara's own copolymers, then predict the LCST of a copolymer nobody reported.

Why it's worth doing:
- Dry lab, laptop only.
- It's squarely my alvaDesc workflow: comonomer descriptors plus a composition number, like the Taguchi project.
- LCST is a well-studied property with a lot of published numbers, so this probably has the most data of all my projects.
- It connects straight to the design question in Ebara's work (what switch temperature will this copolymer have?), and leads naturally into making and testing one on-site (Section 12).

Honest expectation: LCST numbers in the literature are measured different ways (cloud point at 50% transmittance vs onset, different concentrations, heating rates, salt). They're noisy. Molecular weight and end groups also matter. A model that predicts LCST within a few degrees and gets the trends right is the target.

---

## 2. Scope

In scope:
- NIPAM-based copolymers in water with a single comonomer (to start).
- Features: comonomer descriptors (alvaDesc/RDKit), comonomer mole fraction, plus molecular weight and concentration where reported.
- Regression for LCST, trend reproduction, Ebara's systems as a check, one new prediction.

Out of scope:
- Wet lab during the project (Section 12).
- Other thermoresponsive families (PNVCL, POEGMA, etc.) at first. Stretch goal.
- Terpolymers and block copolymers at first. Random copolymers only.
- Predicting drug release or fiber behaviour. Just the LCST.
- Deep learning.

---

## 3. Success criteria

- [ ] A dataset CSV of NIPAM copolymer LCSTs, each row traceable, measurement method recorded.
- [ ] Comonomer descriptors computed and committed as raw files.
- [ ] A regression model with proper validation that beats simple baselines (see 5.4), or a plain statement that it doesn't.
- [ ] Validation that is honest about leakage: tested on comonomers the model never saw, not just new compositions of known comonomers.
- [ ] The model reproduces the direction of the comonomer-hydrophilicity trend and the composition dependence.
- [ ] Predicted LCST for Ebara's copolymer(s) compared to his reported values.
- [ ] One untested copolymer predicted, with an applicability-domain call.
- [ ] A ~2-page note and README in my own voice; one command re-runs it all.

---

## 4. Data specification

### 4.1 Where the data comes from

- Ebara group papers on temperature-responsive nanofibers and smart polymers (SAMURAI publication list, NIMS). These are the check set, not the training set.
- The large literature on NIPAM copolymer LCSTs: review articles on PNIPAM and its copolymers usually tabulate LCST vs comonomer and composition. These are the fastest source.
- Individual papers for comonomer series (acrylamide, acrylic acid at fixed pH, hydroxyethyl methacrylate, butyl methacrylate, etc.).
- Any published cloud-point dataset used in earlier ML papers on LCST, if its data is openly available. Cite it if I use it.

### 4.2 Honest warnings about this data

- Measurement method matters. Cloud point at 50% transmittance, onset of turbidity, and DSC peak can differ by a few degrees. Record the method; model the biggest consistent subset first.
- Concentration, heating rate, and salt shift LCST. Record them.
- Charged comonomers (like acrylic acid) make LCST depend strongly on pH. Either fix the pH or leave those out at first.
- Molecular weight and end groups matter, especially for short chains. Record Mw when reported.
- The biggest trap: data leakage. If the same comonomer appears at five compositions, a random train/test split lets the model "see" that comonomer during training. The honest test is holding out whole comonomers (group k-fold by comonomer). This connects to the leakage-audit work I've been doing.

### 4.3 Dataset schema

- `entry_id`, e.g. `L001`
- `source`, `source_location`, `read_method` (`table` / `digitized_from_figure`)
- `base_monomer` (`NIPAM`)
- `comonomer_name`, `comonomer_smiles` (the comonomer as a monomer molecule, one consistent rule)
- `comonomer_mol_frac` (in the polymer if reported; feed ratio otherwise, and say which in `frac_type`)
- `frac_type` (`polymer_composition` / `feed_ratio`)
- `Mw`, `Mw_method` (blank if not reported)
- `polymer_conc`, `conc_units`
- `solvent` (water or buffer), `pH`, `salt`
- `lcst_method` (`cloud_point_50T` / `onset` / `DSC` / other)
- `heating_rate`
- `lcst_C`
- `notes`

Plus `data/sources.md`.

### 4.4 Expected size

- Likely 80-300 rows across maybe 15-40 different comonomers. The number that matters for the honest test is the number of distinct comonomers, not rows.

Done when:
- [ ] `data/raw/lcst.csv` exists with the schema.
- [ ] I can state rows, distinct comonomers, and rows per `lcst_method`.

---

## 5. Method / pipeline

### 5.1 Stage 1: comonomer structures

- Write SMILES for each comonomer (as the monomer molecule). One rule, written down.
- Check they parse in RDKit, draw them all, look at them.

Done when:
- [ ] `data/processed/comonomers.csv`, all parse, structure grid checked.

### 5.2 Stage 2: descriptors

- alvaDesc on the comonomers during the trial month; commit the raw export.
- RDKit backup set: logP, TPSA, H-bond donors/acceptors, MW, rotatable bonds.
- Join to the dataset. Features per row = comonomer descriptors + mole fraction (+ Mw, concentration where available).
- Also make a simple physically motivated feature: mole fraction x comonomer logP (a rough "how much hydrophobicity did I add"). See if it beats the full descriptor set. It might.

Done when:
- [ ] Descriptor files committed, `model_table.csv` built, no rows lost.

### 5.3 Stage 3: clean features

- Drop constant/near-constant, handle missing, drop one of each highly correlated pair.
- Keep feature count modest relative to the number of distinct comonomers (not rows).

Done when:
- [ ] `features_clean.csv` and a count of surviving features with reasons.

### 5.4 Stage 4: models

Baselines (every real model must beat these):
- Predict 32 C for everything (the PNIPAM value).
- Linear in mole fraction only (ignores which comonomer).

Models:
- Ridge regression on descriptors + composition.
- Random forest / gradient boosting.

Validation, two ways, both reported:
- Random k-fold (optimistic).
- Group k-fold holding out whole comonomers (honest).
The gap between the two is itself a result: it tells me how much the model is just memorising comonomers.

Metrics: MAE in C (easy to interpret), RMSE, R2. Seeds set. No tuning on the reported split.

Done when:
- [ ] `results/metrics.csv` with baselines + models x both validation schemes.
- [ ] A plain sentence about the random vs grouped gap.

### 5.5 Stage 5: reproduce the trends

- Trend 1: at fixed mole fraction, does predicted LCST rise with comonomer hydrophilicity? Rank correlation between predicted LCST and comonomer logP.
- Trend 2: for a comonomer with a published composition series, does predicted LCST follow the measured curve as mole fraction increases?
- Check which features the model leans on, and whether that matches the hydrophilic/hydrophobic explanation.

Done when:
- [ ] Both trends checked with numbers and figures.
- [ ] A paragraph: literature says X; my model says X / doesn't.

### 5.6 Stage 6: Ebara check + new prediction

- Predict LCST for Ebara's reported copolymer(s) at his compositions. Compare to his measured values. Make sure his systems were NOT in the training data for this test (if they were, remove them and refit).
- Predict one untested copolymer: a comonomer/composition pair not in the data. Applicability domain: is the comonomer's descriptor profile inside the training range? If not, say it's an extrapolation.

Done when:
- [ ] Ebara comparison table (predicted vs measured).
- [ ] One new prediction with uncertainty and a domain call.

### 5.7 Stage 7: write it up

- README and a ~2-page note: question, data and its mess, baselines, the random vs grouped result, trends, the Ebara check, the new prediction, limitations.
- My own voice (Section 9).

---

## 6. Repo layout

```
ebara-lcst-model/
  SPEC.md
  README.md
  requirements.txt
  run.sh
  data/
    raw/            <- lcst.csv, alvaDesc export
    processed/
    sources.md
  src/
    s1_comonomers.py
    s2_descriptors.py
    s3_clean.py
    s4_model.py
    s5_trends.py
    s6_ebara_predict.py
  results/
    metrics.csv
    figures/
  notes/
    research_note.md
```

## 7. Tools

- Python 3.11, fresh virtual environment.
- pandas, numpy, scikit-learn (including GroupKFold), matplotlib, rdkit.
- alvaDesc (trial month); commit raw export. RDKit as backup.
- WebPlotDigitizer.
- Git, commit per stage, pinned versions.

## 8. Timeline (rough, part-time)

- Week 1: dataset. Biggest dataset of all my projects, so the most data entry. Start from review tables.
- Week 2: descriptors (alvaDesc trial) and cleaning.
- Week 3: models, both validation schemes, trends.
- Week 4: Ebara check, prediction, write-up.

Note: if all four projects want the alvaDesc month, compute descriptors for all of them in one go. The descriptor step is quick once the SMILES lists are ready.

---

## 9. Writing rules (read this)

Everything a human reads is written by me in my own words: README, the note, code comments, commit messages, emails. Not by an AI, not pasted from one.

Why: the point is to show Ebara (a real person) that I can do this and understand it. AI-polished writing is easy to spot and makes people trust the work less. My own plain writing is more convincing.

How:
- Read, close the source, write from memory.
- Plain and specific. "Adding 10% acrylamide pushed the LCST up by about 4 degrees" beats a vague sentence about tunability.
- Say what I did and what happened, including what didn't work.

Avoid AI tells:
- Words like delve, leverage, robust, comprehensive, seamless, underscore, moreover, furthermore, realm, testament, tapestry, intricate.
- "It's not just X, it's Y" and "it's important to note that."
- Three tidy balanced sentences per paragraph.
- Dashes everywhere, headers on everything.
- Confident filler with no fact.

AI as a tutor, not a ghostwriter:
- Fine: asking it to explain why hydrophobic comonomers lower LCST, or to debug GroupKFold.
- Not fine: pasting its text into anything I hand in or send.
- If AI helps me understand something, I rewrite it myself and make sure I could defend it.
- If AI helps with code, I read and can explain every line, or it doesn't go in.

Before sending anything, read it out loud. If it doesn't sound like me, rewrite it.

---

## 10. Risks

- Mixed measurement methods. Record method, model the biggest consistent subset, add method as a feature only as a secondary analysis.
- Leakage through repeated comonomers. Group k-fold by comonomer is the headline number, not random k-fold.
- Charged comonomers and pH. Fix pH or exclude at first.
- Feed ratio vs actual composition. Record which; prefer actual composition.
- Ebara's systems use comonomers that are rare in the data. Then the Ebara check is an extrapolation; say so rather than hide it.
- Model barely beats the linear baseline. That's a real finding about how much the descriptors add. Report it.

## 11. Stretch goals

- Add a second family (e.g. PNVCL or POEGMA) and see if one model covers both.
- Add salt concentration as a feature for the subset where it varies.
- A small "which comonomer and fraction gives LCST = 38-40 C" design table, the kind of thing useful for a hyperthermia fiber.

---

## 12. The wet half (for the pitch, not for this project)

The natural on-site extension is making the copolymer I predicted (free-radical copolymerisation of NIPAM with the chosen comonomer), measuring its cloud point with a UV-vis spectrometer and a temperature ramp, and seeing how close my prediction was. If there's time, electrospinning it into fibers the way his group does. In the email I should be honest about which of these (if any) I've done in coursework, and say I'd like to learn the rest in his lab. No claims I can't back up.

---

## Changelog

- 2026-09-24: v1 draft written.
- 2026-09-24: Python 3.13 venv (3.11 not installed). Added `role` column to lcst.csv (train / ebara_check) so Ebara rows live in the same file but never enter training.
- 2026-09-26: alvaDesc export (4215 descriptors, 23 molecules) committed as data/raw/alvadesc_export.txt; joined per comonomer in Stage 3 rather than into model_table. Stage 3 filters run on the 22 unique comonomers (not rows). Stage 4: headline subset excludes charged comonomers (AAc, DMAPAAm, DMAPMAm), "all" reported alongside; alvaDesc models pick top-10 fraction-weighted descriptors inside each training fold; grouped CV = GroupKFold(5, shuffled, 5 repeats) plus leave-one-comonomer-out.
