# CRISPRi Flow Cytometry Analysis — 20260928 HLC Trial3 D27 (fixed)

**Revision note:** this version adds compensation (least-squares spillover
matrix from single-stain control wells D1–D4, including a compensated
APC-A channel via the Albumin-AF647 single stain, which reads on the APC-A
detector on this instrument), a cell gate (FSC-A floor + SSC-A upper cap,
a rectangle rather than a floor-only gate), an FSC-H:FSC-A singlet gate, a
time/acquisition-stability QC check, MFI reporting, a per-well
gating-hierarchy table, density-contour scatter plots in place of
overplotted dots, one fixed biexponential cofactor per channel, numbered
plots with Why/How/How-to-read/Result text throughout, and reader-friendly
condition/replicate names in every chart and bar label instead of well IDs.
These changes shift every arm's knockdown numbers from the previous version
of this notebook (see the updated table below); well IDs are preserved only
in metadata/gating-hierarchy tables where provenance matters, and histogram/
scatter traces now pool replicate wells within each condition group rather
than drawing one trace per well.

Finalized marimo notebook covering 6 independent CD81 knockdown arms plus a hepatocyte
(ASGR1/Albumin) marker-stratification analysis, from the fixed HLC Trial3 D27 sample sheet.

## Contents

- `notebook.py` — marimo notebook (open with `marimo edit --sandbox notebook.py`, or
  `uvx marimo edit --sandbox notebook.py`). The PEP 723 header pins dependencies
  (numpy, pandas, plotly, requests, xlrd), so `--sandbox` is required.
- `CONVENTIONS.md` — the chart-design and notebook-organization conventions this analysis
  follows (shared across this line of CRISPRi flow-cytometry notebooks). Read this before
  adding a new arm or chart, so new work matches the existing style.

## Data source

Raw FCS files (one per well) and an `.xls` sample sheet from:
`gs://20240617-landerlab-crispri-project/Flow_Data_CRISPRi/20260928_HLC_Trial3_D27_fixed`

The notebook pulls data live from this GCS bucket via a short-lived OAuth access token
(`gcloud auth print-access-token`, pasted into a password-style widget) — no credentials
are stored in the notebook itself. FCS files are parsed directly (no FlowJo CSV export).
Each well passes through, in order: compensation (inverse spillover matrix built from
single-stain wells D1–D4 against the shared unstained reference A3), a cell gate (FSC-A
floor + SSC-A upper cap, both live sliders), a doublet gate (FSC-Width ceiling, live
slider), and a singlet gate (FSC-H:FSC-A ratio band, the central ~95% of that ratio
computed empirically across all loaded wells after the cell/doublet gates).

## Knockdown metric

CD81 knockdown is reported as a percentage-point shift in "% of gated cells below a fixed
threshold" between a non-targeting guide (ORK/mock) and the CD81-targeting guide, within
each arm's own infection/guide-delivery-gated population. Two complementary gates are
shown side by side throughout: the **1st percentile of the mock population** (captures a
tail/low-signal effect) and the **50th percentile of the mock population** (captures a
broader population-level shift) — arms where the 50th-percentile number is much larger
than the 1st-percentile number indicate a population-wide shift rather than a tail-only
effect. Median fluorescence intensity (MFI) and the MFI ratio vs. the control group are
reported alongside these two percentage-based metrics in every arm's summary table. All
CD81-channel histograms/scatters use a biexponential (`arcsinh(x/cofactor)`) x-axis, with
one fixed cofactor per channel (derived from that channel's signal in the shared unstained
reference well, cached and reused across every plot of that channel) and a left-axis crop
at the 1st percentile of that same unstained reference's own signal.

## The 6 knockdown arms

| Arm | Cell line | Construct | Readout channel | Knockdown, 1st pctile (pp) | Knockdown, 50th pctile (pp) |
|---|---|---|---|---|---|
| 1. pDRT103/pDRT106 split-GFP | 8_3 | All-in-one guide+effector virus requiring `pDRT106` co-infection for GFP reconstitution (GFP+ = confirmed dual infection). Includes a third group, pDRT109 (transduction-only control), not shown in this headline comparison. | BV421-A | 3.9 | 17.6 |
| 2. AA173+AA239 two-virus | 17_3 | Guide-only (`AA173 BFP`) + effector-only (`AA239 dCas9-KRAB-Thy1.1`) viruses, double-positive BFP+/Thy1.1+ gate (Thy1.1 detected via a Thy1.1-FITC antibody). | APC-A | 1.6 | 20.0 |
| 3. AA173+AA239 two-virus | 8_3 | Same construct as arm 2, different cell line. Strongest 1st-percentile knockdown signal of the six arms. | APC-A | 10.4 | 21.6 |
| 4. AA228 mCherry all-in-one | 8_3 | Single-marker (mCherry+) infection gate. Replicates G6/G7 disagree (G6 near-zero/wrong-direction, G7 shows a real effect); the mCherry gate sits on a shoulder, not a clean bimodal valley, a plausible contributor. | BV421-A | 0.3 | -7.3 |
| 5. WTC11-KRAB integrated + AA173 guide | WTC11 | Dox-inducible integrated effector + lentiviral guide; clean 2x2 (guide x Dox). Headline number is Dox-on vs. Dox-off, both within the CD81 guide. | APC-A | 5.8 | 14.9 |
| 6. WTC11 fully-integrated CD81-BFP guide+effector | WTC11 | Effector-delivery comparison: no-effector (BFP+ only) vs. effector-delivered (BFP+ AND Thy1.1+ via AA239), not a Dox comparison. | APC-A | 7.2 | 8.0 |

Each arm's tab includes: a metadata table, a per-well gating-hierarchy table, a pooled
(per condition group) infection-gate-calibrated CD81 histogram + paired density-contour
scatter, reader-friendly per-replicate bar charts (1st- and 50th-percentile-of-mock), and
a tidied summary table with MFI.

## Hepatocyte marker stratification (ASGR1 / Albumin)

Arms 1, 2, 3, 4, and 6 additionally get a within-tab breakdown of CD81 knockdown by
ASGR1+, ALB+, and (arm 1 only) ASGR1+/ALB+ double-positive populations, computed **on
top of** each arm's existing infection/guide-delivery gate (not in place of it). Each
arm shows a top-of-tab strata-comparison bar chart + a data-driven TL;DR, then nested
"All cells"/"ASGR1+"/"ALB+"/"ASGR1+/ALB+" sub-tabs with the full histogram+scatter+bar+
table detail per stratum, all built on the same pooled-trace and reader-friendly-label
conventions as the base arms.

- **ASGR1-FITC is excluded throughout** — that specific antibody clone was independently
  determined to be non-specific/unreliable. Only the alternate **ASGR1-PE** channel
  counts. Arms 1, 2, 3, and 6 have ASGR1-PE available; arm 4 and arm 5 only have the
  excluded ASGR1-FITC clone.
- **Albumin** (goat anti-Albumin + anti-goat Alexa Fluor 647, read in APC-A) is only
  stained in arms 1 and 4.
- **Arm 5 has neither a valid ASGR1 channel nor an Albumin stain** and is excluded from
  the marker-stratified analysis entirely (unstratified only).
- No well in this sample sheet stains the full panel minus ASGR1 or minus Albumin
  specifically, so there is no true FMO control for either marker; ASGR1 gating falls
  back to the fully-unstained per-cell-line naive well, and Albumin gating uses well D1
  (8_3, CD81-BV421 stained but no Albumin) as the closest available partial-stain
  reference — noted as a limitation, not fabricated data.

A separate "Hepatocyte marker gating" calibration section (same live-slider convention
as the main infection-gate calibration section) lets the ASGR1/Albumin thresholds be
adjusted interactively, reactively recomputing every arm's stratified numbers.

## Compensation

Four channels are compensated: BV421-A, FITC-A, PE-A, and APC-A, using single-stain
wells D1 (CD81-BV421), D2 (ASGR1-FITC), D3 (ASGR1-PE), and D4 (Albumin-AF647) against
the shared unstained reference A3 (all four wells are the same 8_3 fixation batch, so
spillover is estimated once and applied notebook-wide, independent of which arm/cell
line a given sample well belongs to). Albumin's Alexa Fluor 647 secondary reads almost
entirely in the APC-A channel on this instrument (confirmed empirically against an
unstained reference), so D4 serves as the APC-A single-stain reference even though no
antibody in this panel is nominally "APC-only." Spillover between these four channels
is small (all off-diagonal coefficients ≤ 0.02 of the primary signal).

## Other methodology notes baked into the analysis

- The Thy1.1-FITC gate (arms 2, 3, 6) uses an **FMO-style stained-but-uninfected**
  reference well per cell line (rather than a fully-unstained naive well), with both
  reference traces shown for comparison in the calibration panel.
- All infection/transduction marker gates (GFP, BFP, Thy1.1, mCherry) are live
  `mo.ui.slider` controls in a dedicated calibration section, and are the actual source
  of truth feeding every arm's knockdown computation — moving a slider reactively
  recomputes the affected arm(s).
- A time/acquisition-stability QC check splits each well's Time channel into 10 bins and
  flags any well whose median FSC-A swings more than 25% between bins (sign of a clog or
  bubble); no well in this dataset was flagged.
- A dedicated investigation found no evidence that the naive-well selection or debris/
  doublet gating contaminated any marker's background reference; apparent
  "bimodal"-looking naive traces turned out to be an artifact of overlaid multi-trace
  plots, not real bimodality in the underlying single-well data (confirmed numerically).
