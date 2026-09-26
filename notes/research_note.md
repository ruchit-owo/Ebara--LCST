<!-- DRAFT written by AI as a starting point. Rewrite in my own voice before anyone reads it (see SPEC.md Section 9). Read it out loud; change anything that doesn't sound like me. -->

# Can a simple model predict where a NIPAM copolymer switches?

## The question

Ebara's temperature-responsive nanofibers are made from NIPAM copolymerised with N-hydroxymethylacrylamide (HMAAm). The HMAAm does two jobs. Its hydroxyl groups crosslink the fibers when they are heated, and because it is more hydrophilic than NIPAM it raises the LCST from about 32 C to wherever the application needs it. For magnetic hyperthermia that is around 42 to 45 C. Choosing the comonomer and how much of it to use is the design decision.

So I asked: given a comonomer's structure and its mole fraction in the chain, how well can the LCST be predicted, and does that still hold for a comonomer the model has never seen?

## Data

I built one table of published LCSTs for random NIPAM copolymers in water or buffer (`data/raw/lcst.csv`). It ended up with 113 training rows from 22 papers, covering 22 comonomers from strongly hydrophilic (acrylamide, HMAAm, hydroxyethylacrylamide) to strongly hydrophobic (butyl methacrylate, styrene, aromatic amino-acid acrylamides). Measured LCSTs run from 15 to 80 C.

The data is messier than the row count suggests.

- LCST was measured in at least five ways. Some papers take the temperature at 50% transmittance, some take the first sign of turbidity, some use an absorbance threshold or a DSC peak. These can differ by a few degrees for the same polymer. I recorded the method for every row but did not correct for it.
- Some papers report the monomer feed and some report the measured composition. For free-radical copolymers these are usually close, but not always.
- 18 of the 22 comonomers come from a single paper each, and 19 rows are values another group compiled from the literature without saying how they were measured.
- Charged comonomers make LCST depend on pH. I left acrylic acid and the two amine comonomers out of the headline analysis.

After removing those, the main set had 103 rows and 19 comonomers.

Ebara's HMAAm copolymers (seven values from four papers) went into a separate check set and were never used for training.

## Method

For each comonomer I calculated RDKit descriptors, and I also ran it through alvaDesc, which gives 4215 descriptors. The feature I thought would matter most was simple: mole fraction times the difference between the comonomer's logP and NIPAM's. Adding 10% of a comonomer that is much more hydrophobic than NIPAM should pull the LCST down, and one that is more hydrophilic should push it up. I did the same with polar surface area and H-bond donors.

Every model had to beat two baselines: 32 C for everything, and a straight line in mole fraction that ignores which comonomer it is.

The part I cared about most was validation. If one comonomer appears at five compositions and the split is random, the model sees that comonomer during training and only has to interpolate. That is not the question. So I report random 5-fold cross-validation next to grouped 5-fold, where all rows of a comonomer are held out together, and leave-one-comonomer-out. The grouped number is the one I trust.

## Results

### Random vs grouped

| Model | Random 5-fold MAE (C) | Grouped 5-fold MAE (C) |
|---|---|---|
| 32 C for everything | 8.9 | 8.9 |
| Linear in mole fraction | 6.2 | 6.6 |
| Ridge, fraction x (logP, TPSA, HBD difference) | 3.6 | 4.3 |
| Random forest, 8 RDKit descriptors | 4.0 | 5.0 |
| Ridge, all 8 RDKit descriptors | 3.8 | 6.8 |
| Ridge, top 10 alvaDesc descriptors (chosen inside each fold) | 4.5 | 20.8 |

Under random splits almost every model looks good, somewhere between 3.6 and 4.5 C. Under grouped splits they separate. The small ridge on four fraction-weighted features only lost 0.7 C. The alvaDesc ridge went from 4.5 to 20.8 C. It had learned the comonomers it was shown and made wild predictions for new ones, mostly for comonomers at high fraction whose descriptor values sat outside what it had seen. I had expected the large descriptor set to help, and it didn't.

Fraction times logP on its own only beat the mole-fraction line by about half a degree. Adding polar surface area and H-bond donors brought the error down from 6.0 to 4.3 C.

![Held-out predictions](../results/figures/parity_logo.png)

### Does it reproduce the known trends?

The literature says hydrophobic comonomers lower the LCST and hydrophilic ones raise it. The data I collected agrees: for each comonomer I fitted the LCST shift per 10 mol%, and it falls with logP (Spearman rho = -0.82 over 19 comonomers). The model, predicting each comonomer at 10 mol% while that comonomer was held out, shows the same order (rho = -0.92).

![logP trend](../results/figures/trend1_logp.png)

For composition, I compared held-out predictions against the measured series. For DMAAm (17 rows from 5 papers) the predicted shift was 5.1 C per 10 mol%, the same as measured. For N-ethylacrylamide it was 3.6 against 3.4. The direction was right for all six series I looked at, but the size was too small for N-tert-butylacrylamide and HIPAAm, and too large for HMAAm (8.1 against 5.1 C per 10 mol%).

![Composition series](../results/figures/trend2_series.png)

One result doesn't fit the textbook picture. The ridge gives polar surface area a negative coefficient, which read on its own says more polar surface lowers the LCST. It is stable (the same sign in 99% of bootstrap refits where I resampled comonomers) and it isn't caused by the bulky aromatic comonomers. My reading is that it is a correction on top of logP: the two fraction-weighted features are correlated at r = -0.74, so the model splits one hydrophilicity effect between them. Taking TPSA out makes grouped error worse (4.3 to 7.2 C), so it is doing something useful, but I wouldn't interpret the coefficients one by one.

### The Ebara check

I removed every HMAAm row from training, so the model had never seen HMAAm, and predicted the seven Ebara group values.

![Ebara check](../results/figures/ebara_check.png)

The mean error was 2.6 C, against 5.8 C for the mole-fraction line and 10.1 C for predicting 32 C. All seven measured values fell inside the 90% prediction interval. The 30 mol% copolymer in water, measured at 56.0 C, was predicted at 56.9 C. Near 20 mol% the model was 3 to 4 C high, mostly for the PBS measurements.

Part of that gap reflects spread in the published HMAAm data rather than the model. Across the papers I used, the LCST rise per 10 mol% HMAAm ranges from about 3 to 7 C, which is normal given the different measurement methods, concentrations and molecular weights. Putting the other groups' HMAAm data back into training only moved the error from 2.6 to 2.1 C.

### A new prediction

I asked which comonomer, not in the data, would put the switch at 38 to 40 C, just above body temperature.

The model predicts that N-acryloylmorpholine (NAM) does this at 14 to 17 mol% (water, 0.1 wt%). At 20 mol% it predicts 41 C, with a 90% interval of 31 to 50 C.

![NAM prediction](../results/figures/new_prediction.png)

NAM's logP, polar surface area and H-bond donor count all lie inside the range of the training comonomers, and its leverage is 0.29 against a threshold of 0.63, so by the descriptors it is inside the model's domain. Structurally it is new: its nearest training comonomer has a fingerprint similarity of only 0.31. I'd call it a prediction inside the descriptor range but for a chemistry the data hasn't covered.

My first choice had been N-(2-hydroxypropyl)methacrylamide (HPMA), because it is hydrophilic and has a hydroxyl group like HMAAm. The model predicts almost no change in LCST for HPMA. Looking at the terms, its logP contribution (+10.7 C at 20 mol%) is cancelled by its TPSA contribution (-9.7 C), which is exactly the correlated term I don't trust. For NAM the prediction comes almost entirely from logP, so I chose it instead.

## Limitations

- 22 comonomers is a small set, and the headline error of 4.3 C comes from 19 of them.
- The mixed measurement methods add noise of a few degrees that no descriptor can explain.
- Molecular weight and end groups are known to matter, especially for short chains, but molecular weight was missing for about half the rows and I couldn't use it.
- Two-dimensional descriptors can't separate isomers. N-n-propylacrylamide gets the same logP as NIPAM even though its homopolymer switches about 10 C lower.
- The model is linear in mole fraction. That works up to about 30 mol% for most comonomers, but the measured curves bend at higher fractions.

## What I'd do next

The obvious test is to make NIPAM-co-NAM at around 15 mol%, measure its cloud point by UV-vis with a temperature ramp, and see whether it lands at 38 to 40 C. I'd also like to learn how HMAAm content and crosslinking change the transition once the polymer is electrospun into fibers, since solution LCST is only the first step for Ebara's systems.
