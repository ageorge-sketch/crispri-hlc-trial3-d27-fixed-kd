# CRISPRi Flow Cytometry Analysis — 20260928 HLC Trial3 D27 (fixed)

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
are stored in the notebook itself. FCS files are parsed directly (no FlowJo CSV export);
debris is gated on FSC-A alone (SSC-A floored at >0, no separate debris criterion), and
doublets/clumps are gated on FSC-Width, both as live adjustable sliders.

## Knockdown metric

CD81 knockdown is reported as a percentage-point shift in "% of gated cells below a fixed
threshold" between a non-targeting guide (ORK/mock) and the CD81-targeting guide, within
each arm's own infection/guide-delivery-gated population. Two complementary gates are
shown side by side throughout: the **1st percentile of the mock population** (captures a
tail/low-signal effect) and the **50th percentile of the mock population** (captures a
broader population-level shift) — arms where the 50th-percentile number is much larger
than the 1st-percentile number indicate a population-wide shift rather than a tail-only
effect. All CD81-channel histograms/scatters use a biexponential (`arcsinh(x/cofactor)`)
x-axis, with a per-arm/per-channel cofactor, so the low-signal region around the gates is
not compressed the way a purely linear or log axis would compress it.

## The 6 knockdown arms

1. **pDRT103/pDRT106 split-GFP, 8_3 cells** — all-in-one guide+effector virus
   (`pDRT103 ORK|CD81`) requiring `pDRT106` co-infection for GFP reconstitution
   (GFP+ = confirmed dual infection). Includes a third group, **pDRT109**
   (transduction-only control, no guide/no effector), which comes out mock-like as
   expected (~0pp on the 1st-percentile metric).
2. **AA173+AA239 two-virus, 17_3 cells** — guide-only (`AA173 BFP`) + effector-only
   (`AA239 dCas9-KRAB-Thy1.1`) viruses, double-positive BFP+/Thy1.1+ gate.
3. **AA173+AA239 two-virus, 8_3 cells** — same construct as arm 2, different cell line,
   kept as its own panel. Strongest 1st-percentile knockdown signal of the six arms.
4. **AA228 mCherry all-in-one, 8_3 cells** — single-marker (mCherry+) infection gate.
   Replicates G6/G7 disagree (G6 near-zero, G7 shows a real effect); the mCherry gate
   sits on a shoulder, not a clean bimodal valley, which is a plausible contributor.
5. **WTC11-KRAB integrated line + AA173 guide, Dox 2x2** — Dox-inducible integrated
   effector + lentiviral guide; clean 2x2 (guide x Dox). Dox-on clearly increases
   knockdown within the CD81 guide, as expected.
6. **WTC11 fully-integrated CD81-BFP guide+effector, C1-C4 only** — redefined mid-session
   from an original Dox=Yes/No design (confounded by an extra virus in the Dox=No wells)
   to a cleaner effector-delivery comparison: no-effector (C3/C4, BFP+ only) vs.
   effector-delivered (C1/C2, BFP+ AND Thy1.1+ via AA239).

Each arm's tab includes: a metadata table, infection-gate-calibrated CD81 histogram +
paired FlowJo-style scatter (one clickable trace per replicate well, never pooled), two
per-replicate bar charts (1st- and 50th-percentile-of-mock), and a tidied summary table.

## Hepatocyte marker stratification (ASGR1 / Albumin)

Arms 1, 2, 3, 4, and 6 additionally get a within-tab breakdown of CD81 knockdown by
ASGR1+, ALB+, and (arm 1 only) ASGR1+/ALB+ double-positive populations, computed **on
top of** each arm's existing infection/guide-delivery gate (not in place of it). Each
arm shows a top-of-tab strata-comparison bar chart + a data-driven TL;DR, then nested
"All cells"/"ASGR1+"/"ALB+"/"ASGR1+/ALB+" sub-tabs with the full histogram+scatter+bar+
table detail per stratum.

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

## Other methodology notes baked into the analysis

- The Thy1.1-FITC gate (arms 2, 3, 6) uses an **FMO-style stained-but-uninfected**
  reference well per cell line (rather than a fully-unstained naive well) — a
  methodology improvement made mid-session, with both reference traces shown for
  comparison in the calibration panel.
- All infection/transduction marker gates (GFP, BFP, Thy1.1, mCherry) are live
  `mo.ui.slider` controls in a dedicated calibration section, and are the actual source
  of truth feeding every arm's knockdown computation — moving a slider reactively
  recomputes the affected arm(s).
- A dedicated investigation found no evidence that the naive-well selection or debris/
  doublet gating contaminated any marker's background reference; apparent
  "bimodal"-looking naive traces turned out to be an artifact of overlaid multi-trace
  plots, not real bimodality in the underlying single-well data (confirmed numerically).
