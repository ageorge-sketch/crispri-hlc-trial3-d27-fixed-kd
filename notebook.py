# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "numpy==2.5.3",
#     "pandas==3.0.6",
#     "plotly==7.1.0",
#     "requests==2.34.2",
#     "xlrd==2.0.2",
# ]
# ///

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium", auto_download=["html"])


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # CRISPRi HLC Trial3 D27 (fixed) -- CD81 knockdown analysis

    Six independent CD81 knockdown arms from the fixed HLC Trial3 D27 sample sheet
    (`20260928_HLC_Trial3_D27_fixed`), each using a different viral-delivery
    construct and/or cell line. Raw FCS files are parsed directly (no FlowJo CSV
    export); events are compensated (least-squares spillover matrix from
    single-stain control wells D1-D4), then gated through a cell gate (FSC-A
    floor + SSC-A upper cap), a doublet gate (FSC-Width ceiling), and a singlet
    gate (FSC-H:FSC-A ratio band) before any infection/guide marker is applied.
    Knockdown is computed as the percentage-point shift in "% of cells below a
    fixed gate" on the compensated CD81 antibody channel between a non-targeting
    guide (ORK) and the CD81-targeting guide, within the infection/guide-positive
    gated population, reported at both the 1st-percentile-of-mock (tail effect)
    and 50th-percentile-of-mock (population-level shift) thresholds, alongside
    median fluorescence intensity (MFI). Every chart and label throughout uses
    reader-friendly condition/replicate names instead of well IDs (well IDs
    appear only in metadata/gating-hierarchy tables where provenance matters).
    Histogram and scatter traces pool replicate wells within each condition
    group; histogram/scatter axes are cropped on the left to the 1st percentile
    of the unstained control well's own signal. ASGR1 and Albumin readouts and
    the K562/unfixed-folder wells are out of scope for the core knockdown
    analysis (ASGR1/Albumin are used only for the marker-stratification tabs
    within each arm) and the K562 wells are excluded throughout.
    """)
    return


@app.cell(hide_code=True)
def _(
    arm1_summary,
    arm2_summary,
    arm3_summary,
    arm4_summary,
    arm5_summary,
    arm6_summary,
    go,
    mo,
    pd,
    plot_overview,
    plot_result,
):
    # Top-level summary intentionally shows only the 6 base unstratified arms
    # (matches the state before the ASGR1/ALB stratification work). The
    # marker-stratified comparison charts/TL;DR/sub-tabs live within each arm's
    # own tab (see build_arm_with_strat) and are not duplicated here.
    _summary_rows = [arm1_summary, arm2_summary, arm3_summary, arm4_summary, arm5_summary, arm6_summary]
    _summary_df = pd.DataFrame(_summary_rows)
    _summary_df["knockdown_pp_1st"] = _summary_df["knockdown_pp_1st"].round(1)
    _summary_df["knockdown_pp_50th"] = _summary_df["knockdown_pp_50th"].round(1)
    _summary_df = _summary_df.rename(columns={
        "arm": "Arm", "cell_line": "Cell line", "mock_group": "Mock group", "guide_group": "Guide-active group",
        "knockdown_pp_1st": "Knockdown, 1st pctile (pp)", "knockdown_pp_50th": "Knockdown, 50th pctile (pp)",
        "flag": "Caveat",
    })
    top_summary_table = mo.ui.table(_summary_df, selection=None)

    _colors = {"1st": "red", "50th": "purple"}
    top_summary_fig = go.Figure()
    top_summary_fig.add_trace(go.Bar(
        x=[r["arm"] for r in _summary_rows], y=[r["knockdown_pp_1st"] for r in _summary_rows],
        name="1st pctile (tail effect)", marker_color=_colors["1st"],
        text=[f"{r['knockdown_pp_1st']:.1f}pp" for r in _summary_rows], textposition="outside",
    ))
    top_summary_fig.add_trace(go.Bar(
        x=[r["arm"] for r in _summary_rows], y=[r["knockdown_pp_50th"] for r in _summary_rows],
        name="50th pctile (population shift)", marker_color=_colors["50th"],
        text=[f"{r['knockdown_pp_50th']:.1f}pp" for r in _summary_rows], textposition="outside",
    ))
    _all_vals = [r["knockdown_pp_1st"] for r in _summary_rows] + [r["knockdown_pp_50th"] for r in _summary_rows]
    _ymax = max(50, max(v for v in _all_vals if v == v) * 1.25) if any(v == v for v in _all_vals) else 50
    top_summary_fig.update_layout(
        title="Plot 1. Bottom line across all 6 arms: knockdown (pp vs. mock), 1st vs. 50th percentile gate",
        barmode="group", yaxis_title="Knockdown (percentage points)", yaxis_range=[0, _ymax],
        height=420, margin=dict(t=60),
        legend=dict(itemclick="toggle", itemdoubleclick="toggleothers"),
    )

    _best1 = max(_summary_rows, key=lambda r: r["knockdown_pp_1st"] if r["knockdown_pp_1st"] == r["knockdown_pp_1st"] else -1e9)
    _best50 = max(_summary_rows, key=lambda r: r["knockdown_pp_50th"] if r["knockdown_pp_50th"] == r["knockdown_pp_50th"] else -1e9)

    top_summary = mo.vstack([
        plot_overview(
            1, "Bottom line across all 6 arms",
            why="States, for every arm, the headline knockdown number at both the tail-effect and population-shift metrics, so the arms can be compared directly.",
            how="Each arm's primary guide-active/effector-active group's % below gate is compared to its own mock/control group's % below gate, at the 1st- and 50th-percentile-of-mock thresholds.",
            how_to_read="Each arm gets one red bar (1st-percentile metric) and one purple bar (50th-percentile metric), both in percentage points of knockdown vs. that arm's own mock group.",
        ),
        top_summary_fig,
        plot_result(
            f"Largest 1st-percentile knockdown: {_best1['arm']} ({_best1['knockdown_pp_1st']:.1f}pp). "
            f"Largest 50th-percentile knockdown: {_best50['arm']} ({_best50['knockdown_pp_50th']:.1f}pp)."
        ),
        top_summary_table,
    ])
    top_summary

    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### How to read each panel below
    Each accordion entry is one construct/cell-line arm: a short description of
    the delivery biology, well metadata and a per-well gating-hierarchy table, a
    peak-normalized pooled histogram + density-contour scatter of the CD81
    readout channel (gated on infection markers, singlet/cell-gated, and
    compensated for cross-channel spillover; dashed red/purple lines = the
    1st-/50th-percentile-of-mock knockdown gates), reader-friendly per-replicate
    bar charts of "% below gate", MFI, and a tidy summary table. Histogram and
    scatter traces pool replicate wells within each condition group; bar charts
    and the gating-hierarchy table always show one entry per physical well,
    labeled by condition group and replicate number rather than a well ID.

    Every plot is numbered (Plot N) with a Why/How/How-to-read note before it and
    a Result line stating the actual computed numbers after it. Numbers follow
    each plot's position in the notebook's own cell-declaration order, not the
    order a reader happens to click through tabs -- a plot inside a later tab can
    have a lower number than one in an earlier tab if it was declared earlier in
    the notebook's dataflow.
    """)
    return


@app.cell(hide_code=True)
def _(
    arm1_content,
    arm2_content,
    arm3_content,
    arm4_content,
    arm5_content,
    arm6_content,
    mo,
):
    mo.ui.tabs({
        "Arm 1. pDRT103/pDRT106 split-GFP, 8_3 cells": arm1_content,
        "Arm 2. AA173+AA239 two-virus, 17_3 cells": arm2_content,
        "Arm 3. AA173+AA239 two-virus, 8_3 cells": arm3_content,
        "Arm 4. AA228 mCherry all-in-one, 8_3 cells": arm4_content,
        "Arm 5. WTC11-KRAB integrated line + AA173 guide, Dox 2x2": arm5_content,
        "Arm 6. WTC11 fully-integrated CD81-BFP guide+effector, Dox-only (use with caution)": arm6_content,
    })
    return


@app.cell(hide_code=True)
def _(mo, pctile_gate):
    # One calibration slider per distinct (channel, reference-well) pair actually
    # used somewhere in the pipeline as an infection/transduction marker gate.
    # The Thy1.1-FITC marker (arms 2/3/6) uses FMO-style stained-but-uninfected
    # reference wells (Thy1.1-FITC stained, no AA239/other virus) instead of the
    # fully-unstained naive, matched per cell line: A8=17_3 (arm 2), D5=8_3
    # (arm 3), A4=WTC11 (arm 6) -- distinct keys (suffixed _fmo) from the GFP/BFP
    # gates so switching the Thy1.1 reference doesn't also move arm 1's GFP gate
    # or arms 2/3/5/6's BFP gates, even where they'd otherwise share a channel.
    _infection_gate_defaults = {
        "FITC-A__A3": pctile_gate("A3", "FITC-A", 99),
        "BV421-A__A2": pctile_gate("A2", "BV421-A", 99),
        "FITC-A__A8_fmo": pctile_gate("A8", "FITC-A", 99),
        "BV421-A__A3": pctile_gate("A3", "BV421-A", 99),
        "FITC-A__D5_fmo": pctile_gate("D5", "FITC-A", 99),
        "Y610-mCHERRY-A__A3": pctile_gate("A3", "Y610-mCHERRY-A", 99),
        "BV421-A__A1": pctile_gate("A1", "BV421-A", 99),
        "FITC-A__A4_fmo": pctile_gate("A4", "FITC-A", 99),
        # Hepatocyte markers: ASGR1 (non-FITC channel only -- that FITC clone is
        # known non-specific per the user, excluded entirely) and Albumin
        # (Alexa Fluor 647 secondary, read in the APC-A channel -- free of CD81
        # readout conflict since no well stains both CD81-APC and Albumin-AF647).
        # ASGR1-PE naive reference, per cell line (no true FMO well exists for
        # ASGR1/Albumin in this sample sheet -- falling back to fully-unstained).
        "PE-A__A1": pctile_gate("A1", "PE-A", 99),
        "PE-A__A2": pctile_gate("A2", "PE-A", 99),
        "PE-A__A3": pctile_gate("A3", "PE-A", 99),
        # Albumin-AF647/APC-A reference: D1 (8_3, CD81-BV421 stained, no
        # Albumin) used instead of the fully-unstained A3, since it at least
        # matches the BV421-spillover background of the arm-1/arm-4 sample wells
        # -- the closest available approximation to a true FMO for this marker.
        "APC-A__D1_alb": pctile_gate("D1", "APC-A", 99),
    }

    def _widget_label(key):
        channel, ref = key.split("__")
        if ref.endswith("_fmo"):
            return f"{channel} infection-gate threshold (FMO ref = well {ref.replace('_fmo', '')}, Thy1.1-stained/uninfected)"
        if ref.endswith("_alb"):
            base_well = ref.replace("_alb", "")
            return f"{channel} Albumin-gate threshold (partial-stain ref = well {base_well}, CD81-stained/no Albumin -- no true FMO available)"
        return f"{channel} infection-gate threshold (naive/unstained ref = well {ref})"

    infection_gate_widgets = mo.ui.dictionary({
        key: mo.ui.slider(
            start=0, stop=max(2000.0, round(val * 6)), step=max(50.0, round(val / 100)),
            value=round(val), show_value=True, full_width=True,
            label=_widget_label(key),
        )
        for key, val in _infection_gate_defaults.items()
    })
    infection_gate_widgets
    return (infection_gate_widgets,)


@app.cell(hide_code=True)
def _(
    biexp,
    biexp_ticks,
    calibration_trace_labels,
    debris_gate,
    gated_range,
    infection_gate_widgets,
    interactive_hist,
    mo,
    np,
    plot_overview,
    plot_result,
    sample_sheet,
):
    _CALIBRATION_CFGS = [
        {"key": "FITC-A__A3", "channel": "FITC-A", "naive": "A3", "plot_n": 2,
         "used_by": {"Arm 1 GFP marker (B1-B6)": ["B1", "B2", "B3", "B4", "B5", "B6"]}},
        {"key": "BV421-A__A2", "channel": "BV421-A", "naive": "A2", "plot_n": 3,
         "used_by": {"Arm 2 BFP marker (C5-C8)": ["C5", "C6", "C7", "C8"]}},
        {"key": "FITC-A__A8_fmo", "channel": "FITC-A", "naive": "A8", "secondary_naive": "A2", "plot_n": 4,
         "used_by": {"Arm 2 Thy1.1 marker (C5-C8)": ["C5", "C6", "C7", "C8"]}},
        {"key": "BV421-A__A3", "channel": "BV421-A", "naive": "A3", "plot_n": 5,
         "used_by": {"Arm 3 BFP marker (C9-C12)": ["C9", "C10", "C11", "C12"]}},
        {"key": "FITC-A__D5_fmo", "channel": "FITC-A", "naive": "D5", "secondary_naive": "A3", "plot_n": 6,
         "used_by": {"Arm 3 Thy1.1 marker (C9-C12)": ["C9", "C10", "C11", "C12"]}},
        {"key": "Y610-mCHERRY-A__A3", "channel": "Y610-mCHERRY-A", "naive": "A3", "plot_n": 7,
         "used_by": {"Arm 4 mCherry marker (G4-G7)": ["G4", "G5", "G6", "G7"]}},
        {"key": "BV421-A__A1", "channel": "BV421-A", "naive": "A1", "plot_n": 8,
         "used_by": {"Arm 5 BFP marker (H3-H8)": ["H3", "H4", "H5", "H6", "H7", "H8"],
                     "Arm 6 BFP marker (C1-C4)": ["C1", "C2", "C3", "C4"]}},
        {"key": "FITC-A__A4_fmo", "channel": "FITC-A", "naive": "A4", "secondary_naive": "A1", "plot_n": 9,
         "used_by": {"Arm 6 Thy1.1 marker (C1/C2)": ["C1", "C2"]}},
    ]

    def build_calibration_panel(cfg):
        channel, naive, key = cfg["channel"], cfg["naive"], cfg["key"]
        slider = infection_gate_widgets.elements[key]
        gate = slider.value
        naive_vals = debris_gate(naive)[channel].values
        naive_cell_line = sample_sheet.loc[naive, "Cell Line"]
        is_fmo = key.endswith("_fmo")
        is_alb = key.endswith("_alb") or cfg.get("ref_kind") == "alb"
        if is_fmo:
            naive_label = f"FMO ref ({naive_cell_line}, Thy1.1-stained/uninfected, well {naive})"
        elif is_alb:
            naive_label = f"Partial-stain ref ({naive_cell_line}, CD81-stained/no Albumin, well {naive} -- no true FMO available)"
        else:
            naive_label = f"Naive ({naive_cell_line}, unstained)"

        traces = {naive_label: naive_vals}
        pooled_for_range = [naive_vals]

        # Keep the fully-unstained trace visible alongside an FMO reference so the
        # user can see the difference directly, per the FMO-methodology request.
        secondary = cfg.get("secondary_naive")
        if secondary:
            sec_vals = debris_gate(secondary)[channel].values
            sec_cell_line = sample_sheet.loc[secondary, "Cell Line"]
            sec_label = f"(for comparison) Naive ({sec_cell_line}, unstained, well {secondary})"
            traces[sec_label] = sec_vals
            pooled_for_range.append(sec_vals)

        # Per-well traces (not pooled per-arm) labeled with readable
        # construct/condition text -- never a bare well ID in the legend -- via
        # the shared calibration_trace_labels helper.
        all_sample_wells = [w for wells in cfg["used_by"].values() for w in wells]
        well_labels = calibration_trace_labels(all_sample_wells)

        for w in all_sample_wells:
            vals = debris_gate(w)[channel].values
            traces[well_labels[w]] = vals
            pooled_for_range.append(vals)
        pooled = np.concatenate(pooled_for_range)
        xr = gated_range(naive_vals, pooled, hi_pct=99.0)
        cofactor = max(gate, 1.0) / 5.0

        def _tx(v):
            return biexp(v, cofactor)

        traces_t = {k: _tx(v) for k, v in traces.items()}
        xr_t = (float(_tx(xr[0])), float(_tx(xr[1])))
        gate_t = float(_tx(gate))
        _plot_n = cfg.get("plot_n", 1)
        fig = interactive_hist(
            traces_t, gate_t, f"gate = {gate:,.0f}", xr_t,
            f"Plot {_plot_n}. {channel} infection-gate calibration (reference = {naive_label})",
        )
        tickvals, ticktext = biexp_ticks(xr[0], xr[1], cofactor)
        fig.update_layout(xaxis=dict(tickvals=tickvals, ticktext=ticktext, title=f"{channel} (biexponential scale)"))

        _frac_above = float(np.mean(np.concatenate([debris_gate(w)[channel].values for w in all_sample_wells]) > gate)) if all_sample_wells else float("nan")
        blocks = [
            plot_overview(
                _plot_n, f"{channel} infection-gate calibration",
                why=f"Shows where the {channel} infection/transduction-marker gate sits relative to the reference well(s) and the actual sample wells that use it.",
                how="Reference well(s) and every sample well using this (channel, reference) pair are gated (singlet/cell-gated, compensated) and plotted as biexponential histograms; the slider sets the gate threshold live.",
                how_to_read="Each line is one well (reference or sample), labeled by condition, never a bare well ID. The dashed red line is the current gate threshold.",
            ),
            mo.md(f"**Used by:** {', '.join(cfg['used_by'].keys())}"), slider, fig,
            plot_result(f"At the current gate, {_frac_above:.1%} of pooled events across this gate's sample wells fall above the threshold."),
        ]

        if is_fmo:
            blocks.append(mo.md(
                "*FMO-style reference: this gate is now calibrated off a "
                "**stained-but-uninfected** control (Thy1.1-FITC antibody applied, "
                "no AA239/other virus delivered) matched to this arm's own cell "
                "line, rather than the fully-unstained naive well -- a methodology "
                "improvement agreed with the user, since a true FMO control better "
                "isolates antibody/channel background from biology than a "
                "completely unstained well does. The fully-unstained trace is kept "
                "visible above for comparison.*"
            ))

        if key == "Y610-mCHERRY-A__A3":
            _frac = {w: float(np.mean(debris_gate(w)[channel].values > gate)) for w in ["G4", "G5", "G6", "G7"]}
            blocks.append(mo.callout(mo.md(
                "**Arm 4 mCherry gate check (requested):** the mCherry distributions in "
                "G4-G7 show no clean bimodal valley -- percentiles rise smoothly from the "
                "10th to 99th with no plateau/dip, so this 99th-pctile-of-unstained "
                "threshold sits on a **shoulder, not a valley** (only ~10-18% of cells "
                "clear it in any well, i.e. it is not splitting off an obvious separate "
                "positive population). At the current gate, mCherry+ fractions are "
                f"G4={_frac['G4']:.1%}, G5={_frac['G5']:.1%}, G6={_frac['G6']:.1%}, "
                f"G7={_frac['G7']:.1%} -- **G7 has ~1.8x the mCherry+ fraction of G6** "
                "despite both being CD81-guide replicates. This differential "
                "infection/transduction efficiency between G6 and G7, combined with an "
                "ambiguous (non-bimodal) gate placement, is a plausible contributor to "
                "why G6 and G7 disagreed in the knockdown metrics (G6 near-zero/"
                "wrong-direction, G7 showing a real effect) -- though it does not by "
                "itself prove the arm 4 population-level (-6.5pp) result is purely a "
                "gating artifact; treat arm 4 with extra caution either way."
            ), kind="warn"))

        return mo.vstack(blocks)

    infection_calibration_accordion = mo.ui.tabs({
        (
            f"{cfg['channel']} FMO ref: {sample_sheet.loc[cfg['naive'], 'Cell Line']} (Thy1.1-stained/uninfected)"
            if cfg["key"].endswith("_fmo")
            else f"{cfg['channel']} naive: {sample_sheet.loc[cfg['naive'], 'Cell Line']} (unstained)"
        ): build_calibration_panel(cfg)
        for cfg in _CALIBRATION_CFGS
    })

    infection_calibration_section = mo.vstack([
        mo.md(
            "## Infection-gate calibration\n"
            "FlowJo-style gating histograms (biexponential x-axis, same "
            "`arcsinh(x/cofactor)` transform used for the CD81 readout channels above) "
            "for every infection/transduction marker gate actually used in the pipeline. "
            "The Thy1.1-FITC gate (arms 2/3/6) uses an FMO-style stained-but-uninfected "
            "reference well per cell line rather than the fully-unstained naive (both "
            "traces are shown for comparison). **The sliders below are now the live "
            "source of truth** for these gates -- moving one reactively re-runs every "
            "arm that uses that (channel, reference-well) pair and updates its "
            "knockdown numbers above."
        ),
        infection_calibration_accordion,
    ])
    infection_calibration_section
    return (build_calibration_panel,)


@app.cell(hide_code=True)
def _(build_calibration_panel, mo, sample_sheet):
    # Hepatocyte marker gating (ASGR1 non-FITC channel only + Albumin), reusing
    # the same build_calibration_panel/infection_gate_widgets machinery as the
    # infection-gate section above. ASGR1-FITC is excluded everywhere (per the
    # user, that clone is non-specific) -- only ASGR1-PE (read in PE-A) counts.
    _HEP_MARKER_CFGS = [
        {"key": "PE-A__A3", "channel": "PE-A", "naive": "A3", "plot_n": 10,
         "used_by": {"Arm 1 ASGR1 marker (B3-B6)": ["B3", "B4", "B5", "B6"],
                     "Arm 3 ASGR1 marker (C9-C12)": ["C9", "C10", "C11", "C12"]}},
        {"key": "PE-A__A2", "channel": "PE-A", "naive": "A2", "plot_n": 11,
         "used_by": {"Arm 2 ASGR1 marker (C5-C8)": ["C5", "C6", "C7", "C8"]}},
        {"key": "PE-A__A1", "channel": "PE-A", "naive": "A1", "plot_n": 12,
         "used_by": {"Arm 6 ASGR1 marker (C1-C4)": ["C1", "C2", "C3", "C4"]}},
        {"key": "APC-A__D1_alb", "channel": "APC-A", "naive": "D1", "ref_kind": "alb", "plot_n": 13,
         "used_by": {"Arm 1 Albumin marker (B3-B6)": ["B3", "B4", "B5", "B6"],
                     "Arm 4 Albumin marker (G4-G7)": ["G4", "G5", "G6", "G7"]}},
    ]

    hepatocyte_marker_tabs = mo.ui.tabs({
        (
            f"{cfg['channel']} ASGR1 naive: {sample_sheet.loc[cfg['naive'], 'Cell Line']}"
            if cfg["key"].startswith("PE-A")
            else f"{cfg['channel']} Albumin ref: {sample_sheet.loc[cfg['naive'], 'Cell Line']} (partial-stain)"
        ): build_calibration_panel(cfg)
        for cfg in _HEP_MARKER_CFGS
    })

    hepatocyte_marker_section = mo.vstack([
        mo.md(
            "## Hepatocyte marker gating (ASGR1, Albumin)\n"
            "**ASGR1-FITC is excluded throughout this section and the stratified "
            "knockdown analysis below** -- per the user, that specific ASGR1-FITC "
            "antibody/clone is known non-specific/unreliable. Only wells stained with "
            "the alternate **ASGR1-PE** channel are used for ASGR1 gating: arms 1, 2, "
            "3, and 6 have ASGR1-PE available (gated below); **arm 5 (H3-H8) and arm "
            "4 (G4-G7) only have ASGR1-FITC and are excluded from ASGR1 gating** "
            "(arm 4 still gets Albumin gating, arm 5 has neither marker and is "
            "dropped from the stratified section entirely). **Albumin** (goat "
            "anti-Albumin + anti-goat Alexa Fluor 647, read in APC-A) is only "
            "stained in arms 1 and 4. No well in this sample sheet stains the full "
            "panel minus ASGR1 or minus Albumin specifically, so there is no true "
            "FMO control for either marker; ASGR1 falls back to the fully-unstained "
            "per-cell-line naive, and Albumin uses well D1 (8_3, CD81-BV421 stained "
            "but no Albumin) as the closest available partial-stain reference -- "
            "noted as a limitation, not fabricated data."
        ),
        hepatocyte_marker_tabs,
    ])
    hepatocyte_marker_section
    return


@app.cell(hide_code=True)
def _(
    MIN_CELLS_PER_GATE,
    apply_markers,
    infection_gate_widgets,
    is_control_label,
    np,
    pd,
):
    # Marker-stratified CD81 knockdown: same 1st/50th-percentile-of-mock metric
    # as every arm above, but computed on top of (AND'd with) each arm's existing
    # infection/guide-delivery gate -- never in place of it. Arm 5 is excluded
    # entirely (no valid ASGR1 channel, no Albumin stain at all on H3-H8).
    def stratum_filter(d, asgr1_spec, alb_spec):
        if asgr1_spec is not None:
            d = d[d[asgr1_spec[0]].values > asgr1_spec[1]]
        if alb_spec is not None:
            d = d[d[alb_spec[0]].values > alb_spec[1]]
        return d

    def stratified_knockdown_rows(arm_label, groups, readout_channel, control_group_index, asgr1_gate=None, alb_gate=None):
        strata = [("All cells (unstratified)", None, None)]
        if asgr1_gate is not None:
            strata.append(("ASGR1+", ("PE-A", asgr1_gate), None))
        if alb_gate is not None:
            strata.append(("ALB+", None, ("APC-A", alb_gate)))
        if asgr1_gate is not None and alb_gate is not None:
            strata.append(("ASGR1+/ALB+", ("PE-A", asgr1_gate), ("APC-A", alb_gate)))

        control_label, control_wells, _ = groups[control_group_index]
        rows = []
        for stratum_label, asgr1_spec, alb_spec in strata:
            infection_gated = {}
            for label, wells, marker_specs in groups:
                for w in wells:
                    d = apply_markers(w, marker_specs)  # existing infection/guide gate, unchanged
                    d = stratum_filter(d, asgr1_spec, alb_spec)  # AND'd marker stratification
                    infection_gated[w] = (label, d)

            control_vals = np.concatenate([infection_gated[w][1][readout_channel].values for w in control_wells]) \
                if all(len(infection_gated[w][1]) > 0 for w in control_wells) else np.array([])
            gate1 = float(np.percentile(control_vals, 1)) if len(control_vals) else None
            gate50 = float(np.percentile(control_vals, 50)) if len(control_vals) else None

            stratum_rows = []
            for label, wells, _ in groups:
                for w in wells:
                    d = infection_gated[w][1]
                    n = len(d)
                    if n < MIN_CELLS_PER_GATE and not is_control_label(label):
                        stratum_rows.append({"arm": arm_label, "stratum": stratum_label, "well": w, "label": label, "n": n,
                                              "pct_below_1st": np.nan, "pct_below_50th": np.nan, "excluded_low_n": True})
                        continue
                    p1 = 100 * np.mean(d[readout_channel].values < gate1) if gate1 is not None and n else np.nan
                    p50 = 100 * np.mean(d[readout_channel].values < gate50) if gate50 is not None and n else np.nan
                    stratum_rows.append({"arm": arm_label, "stratum": stratum_label, "well": w, "label": label, "n": n,
                                          "pct_below_1st": p1, "pct_below_50th": p50, "excluded_low_n": False})

            _sdf = pd.DataFrame(stratum_rows)
            ctrl_mean_1 = _sdf.loc[_sdf["label"] == control_label, "pct_below_1st"].mean()
            ctrl_mean_50 = _sdf.loc[_sdf["label"] == control_label, "pct_below_50th"].mean()
            _sdf["knockdown_pp_1st"] = _sdf["pct_below_1st"] - ctrl_mean_1
            _sdf["knockdown_pp_50th"] = _sdf["pct_below_50th"] - ctrl_mean_50
            rows.extend(_sdf.to_dict("records"))
        return rows

    STRAT_ARM_CFGS = [
        {
            "arm_label": "Arm 1: pDRT103/106, 8_3", "plot_start": 43, "naive_well": "A3",
            "groups": [
                ("ORK (pDRT103 ORK + pDRT106)", ["B3", "B4"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
                ("pDRT109 (transduction-only control + pDRT106)", ["B1", "B2"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
                ("CD81 (pDRT103 CD81 + pDRT106)", ["B5", "B6"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
            ],
            "readout_channel": "BV421-A", "control_group_index": 0,
            "asgr1_gate": infection_gate_widgets.value["PE-A__A3"],
            "alb_gate": infection_gate_widgets.value["APC-A__D1_alb"],
        },
        {
            "arm_label": "Arm 2: AA173+AA239, 17_3", "plot_start": 56, "naive_well": "A2",
            "groups": [
                ("ORK (AA173 ORK + AA239)", ["C5", "C6"],
                 [("BV421-A", infection_gate_widgets.value["BV421-A__A2"]), ("FITC-A", infection_gate_widgets.value["FITC-A__A8_fmo"])]),
                ("CD81 (AA173 CD81 + AA239)", ["C7", "C8"],
                 [("BV421-A", infection_gate_widgets.value["BV421-A__A2"]), ("FITC-A", infection_gate_widgets.value["FITC-A__A8_fmo"])]),
            ],
            "readout_channel": "APC-A", "control_group_index": 0,
            "asgr1_gate": infection_gate_widgets.value["PE-A__A2"],
            "alb_gate": None,
        },
        {
            "arm_label": "Arm 3: AA173+AA239, 8_3", "plot_start": 61, "naive_well": "A3",
            "groups": [
                ("ORK (AA173 ORK + AA239)", ["C9", "C10"],
                 [("BV421-A", infection_gate_widgets.value["BV421-A__A3"]), ("FITC-A", infection_gate_widgets.value["FITC-A__D5_fmo"])]),
                ("CD81 (AA173 CD81 + AA239)", ["C11", "C12"],
                 [("BV421-A", infection_gate_widgets.value["BV421-A__A3"]), ("FITC-A", infection_gate_widgets.value["FITC-A__D5_fmo"])]),
            ],
            "readout_channel": "APC-A", "control_group_index": 0,
            "asgr1_gate": infection_gate_widgets.value["PE-A__A3"],
            "alb_gate": None,
        },
        {
            "arm_label": "Arm 4: AA228 mCherry, 8_3", "plot_start": 66, "naive_well": "A3",
            "groups": [
                ("ORK (AA228 mCherry ORK)", ["G4", "G5"], [("Y610-mCHERRY-A", infection_gate_widgets.value["Y610-mCHERRY-A__A3"])]),
                ("CD81 (AA228 mCherry CD81)", ["G6", "G7"], [("Y610-mCHERRY-A", infection_gate_widgets.value["Y610-mCHERRY-A__A3"])]),
            ],
            "readout_channel": "BV421-A", "control_group_index": 0,
            "asgr1_gate": None,  # arm 4 only has ASGR1-FITC (excluded)
            "alb_gate": infection_gate_widgets.value["APC-A__D1_alb"],
        },
        {
            "arm_label": "Arm 6: WTC11 fully-integrated, C1-C4", "plot_start": 71, "naive_well": "A1",
            "groups": [
                ("No effector (BFP+ only, C3/C4)", ["C3", "C4"], [("BV421-A", infection_gate_widgets.value["BV421-A__A1"])]),
                ("Effector delivered (BFP+ AND Thy1.1+, C1/C2)", ["C1", "C2"],
                 [("BV421-A", infection_gate_widgets.value["BV421-A__A1"]), ("FITC-A", infection_gate_widgets.value["FITC-A__A4_fmo"])]),
            ],
            "readout_channel": "APC-A", "control_group_index": 0,
            "asgr1_gate": infection_gate_widgets.value["PE-A__A1"],
            "alb_gate": None,
        },
        # Arm 5 (H3-H8) intentionally omitted: ASGR1-FITC only (excluded) and no
        # Albumin stain at all on these wells -- neither marker is available.
    ]

    _strat_all_rows = []
    for _cfg in STRAT_ARM_CFGS:
        _strat_all_rows.extend(stratified_knockdown_rows(
            _cfg["arm_label"], _cfg["groups"], _cfg["readout_channel"], _cfg["control_group_index"],
            asgr1_gate=_cfg["asgr1_gate"], alb_gate=_cfg["alb_gate"],
        ))

    marker_strat_df = pd.DataFrame(_strat_all_rows)
    marker_strat_df
    return STRAT_ARM_CFGS, marker_strat_df, stratum_filter


@app.cell(hide_code=True)
def _(
    apply_markers,
    biexp,
    channel_cofactor,
    debris_gate,
    gated_range,
    go,
    interactive_hist,
    marker_strat_df,
    mo,
    np,
    pd,
    plot_overview,
    plot_result,
    replicate_labels,
    scatter_gate,
    stratum_filter,
):
    def _tidy_strat_table(df):
        _d = df.copy()
        for c in ["pct_below_1st", "pct_below_50th", "knockdown_pp_1st", "knockdown_pp_50th"]:
            _d[c] = _d[c].round(1)
        _d = _d.drop(columns=["excluded_low_n"])
        _d = _d.rename(columns={
            "stratum": "Stratum", "well": "Well", "label": "Condition", "n": "n (gated cells)",
            "pct_below_1st": "% below 1st-pctile-mock", "pct_below_50th": "% below 50th-pctile-mock",
            "knockdown_pp_1st": "Knockdown, 1st pctile (pp)", "knockdown_pp_50th": "Knockdown, 50th pctile (pp)",
        })
        return mo.ui.table(_d, selection=None)

    def _strat_bar_chart(arm_label, stratum_label, metric_col, title_suffix, plot_n, groups):
        _sub = marker_strat_df[(marker_strat_df["arm"] == arm_label) & (marker_strat_df["stratum"] == stratum_label)]
        if not len(_sub):
            return None
        rl = {}
        for label, wells, _ in groups:
            rl.update(replicate_labels([w for w in wells if w in _sub["well"].values], group_label=label))
        colors_map = {g[0]: c for g, c in zip(groups, ["#4C78A8", "#E45756", "#54A24B", "#F58518"])}
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=[rl.get(w, w) for w in _sub["well"]], y=_sub[metric_col],
            text=[f"{v:.1f}%<br>(n={n:,})" for v, n in zip(_sub[metric_col], _sub["n"])],
            textposition="outside",
            marker_color=[colors_map.get(l, "#999") for l in _sub["label"]],
        ))
        ymax = max(50, float(_sub[metric_col].max()) * 1.25) if _sub[metric_col].notna().any() else 100
        fig.update_layout(
            title=f"Plot {plot_n}. {arm_label} -- {stratum_label}: {title_suffix} (per replicate)",
            yaxis_title="% below gate", yaxis_range=[0, min(100, ymax)], height=360, margin=dict(t=60),
        )
        return fig

    def build_marker_stratum_block(cfg, stratum_label, asgr1_spec, alb_spec, plot_start):
        """Content for ONE stratum (e.g. 'ASGR1+') of one arm: pooled biexponential
        CD81 histogram + paired density scatter (same convention as the main arms
        section, dual 1st-/50th-percentile-of-mock dashed threshold lines) +
        reader-friendly per-replicate bar charts + a table filtered to this
        stratum."""
        groups = cfg["groups"]
        readout_channel = cfg["readout_channel"]
        control_group_index = cfg["control_group_index"]
        naive_well = cfg["naive_well"]
        arm_label = cfg["arm_label"]
        primary_label = groups[-1][0]
        control_label, control_wells, _ = groups[control_group_index]
        naive_vals_raw = debris_gate(naive_well)[readout_channel].values

        gated_by_well = {}
        for label, wells, marker_specs in groups:
            for w in wells:
                d = apply_markers(w, marker_specs)
                d = stratum_filter(d, asgr1_spec, alb_spec)
                gated_by_well[w] = (label, d)

        control_vals = np.concatenate([gated_by_well[w][1][readout_channel].values for w in control_wells]) \
            if all(len(gated_by_well[w][1]) > 0 for w in control_wells) else np.array([])
        if not len(control_vals):
            return mo.md(f"**{stratum_label}:** too few control cells to draw a gate -- skipped.")
        gate1 = float(np.percentile(control_vals, 1))
        gate50 = float(np.percentile(control_vals, 50))

        pooled_vals = np.concatenate([d[readout_channel].values for _, d in gated_by_well.values() if len(d)])
        xr = gated_range(naive_vals_raw, pooled_vals, hi_pct=99.0)
        cofactor = channel_cofactor(readout_channel)

        def _tx(v):
            return biexp(v, cofactor)

        xr_t = (float(_tx(xr[0])), float(_tx(xr[1])))
        gate1_t, gate50_t = float(_tx(gate1)), float(_tx(gate50))

        # Histogram and scatter traces pool replicate wells within each group.
        hist_by_group = {}
        scatter_by_group = {}
        for w, (label, d) in gated_by_well.items():
            if len(d):
                vals_t = _tx(d[readout_channel].values)
                hist_by_group.setdefault(label, []).append(vals_t)
                xs, ys = scatter_by_group.setdefault(label, ([], []))
                xs.append(vals_t)
                ys.append(d["SSC-A"].values)
        hist_traces = {"unstained reference": _tx(naive_vals_raw)}
        hist_traces.update({label: np.concatenate(vs) for label, vs in hist_by_group.items()})
        scatter_traces = {
            label: (np.concatenate(xs), np.concatenate(ys)) for label, (xs, ys) in scatter_by_group.items()
        }

        fig_hist = interactive_hist(
            hist_traces, gate1_t, f"1st %ile mock ({control_label}) = {gate1:,.0f}",
            xr_t, f"Plot {plot_start}. {arm_label} -- {stratum_label}: {readout_channel} distribution",
            gate2=gate50_t, gate2_label=f"50th %ile mock ({control_label}) = {gate50:,.0f}", gate2_color="purple",
        )
        fig_scatter = scatter_gate(
            scatter_traces, xr_t, f"Plot {plot_start + 1}. {arm_label} -- {stratum_label}: {readout_channel} vs SSC-A",
            y_chan="SSC-A", gate=gate1_t, gate_label="1st %ile mock", gate2=gate50_t, gate2_label="50th %ile mock", gate2_color="purple",
        )
        bar1 = _strat_bar_chart(arm_label, stratum_label, "pct_below_1st", "% below 1st-pctile-mock gate", plot_start + 2, groups)
        bar50 = _strat_bar_chart(arm_label, stratum_label, "pct_below_50th", "% below 50th-pctile-mock gate", plot_start + 3, groups)

        _strat_sub = marker_strat_df[(marker_strat_df["arm"] == arm_label) & (marker_strat_df["stratum"] == stratum_label)]
        _kd1 = float(_strat_sub.loc[_strat_sub["label"] == primary_label, "knockdown_pp_1st"].mean()) if len(_strat_sub) else float("nan")
        _kd50 = float(_strat_sub.loc[_strat_sub["label"] == primary_label, "knockdown_pp_50th"].mean()) if len(_strat_sub) else float("nan")
        _strat_n = int(_strat_sub.loc[_strat_sub["label"] == primary_label, "n"].sum()) if len(_strat_sub) else 0

        blocks = [
            plot_overview(
                plot_start, f"{arm_label} -- {stratum_label}: {readout_channel} distribution",
                why=f"Checks whether the {primary_label} vs. {control_label} knockdown signal seen across all cells in {arm_label} still holds within the {stratum_label} cell group specifically.",
                how=f"Same computation as the arm's unstratified histogram, restricted to cells additionally passing the {stratum_label} marker gate, pooled by condition group.",
                how_to_read="Same layout as the unstratified histogram: biexponential x-axis, % of each group's peak count on the y-axis, dashed red/purple lines mark the 1st- and 50th-percentile-of-mock gates for this stratum.",
            ),
            fig_hist,
            plot_overview(
                plot_start + 1, f"{arm_label} -- {stratum_label}: {readout_channel} vs SSC-A",
                why=f"Checks the {stratum_label} population separation in 2D, the same purpose as the arm's unstratified scatter plot.",
                how="Same computation as the arm's unstratified scatter plot, restricted to this stratum.",
                how_to_read="Same density-contour layout as the unstratified scatter plot above.",
            ),
            fig_scatter,
        ]
        if bar1 is not None:
            blocks.append(plot_overview(
                plot_start + 2, f"{arm_label} -- {stratum_label}: % below 1st-percentile-of-mock gate",
                why=f"States the tail-effect knockdown number for the {stratum_label} stratum directly, per replicate well.",
                how="Same computation as the arm's unstratified 1st-percentile bar chart, restricted to this stratum.",
                how_to_read="Same layout as the arm's unstratified bar chart; bars are labeled by condition group and replicate number, never a well ID.",
            ))
            blocks.append(bar1)
            blocks.append(plot_result(
                f"Mean knockdown for {primary_label} vs. {control_label} within the {stratum_label} stratum "
                f"is {_kd1:+.1f} percentage points (n={_strat_n:,} {primary_label} cells)."
            ))
        if bar50 is not None:
            blocks.append(plot_overview(
                plot_start + 3, f"{arm_label} -- {stratum_label}: % below 50th-percentile-of-mock gate",
                why=f"States the population-level-shift knockdown number for the {stratum_label} stratum directly, per replicate well.",
                how="Same computation as the arm's unstratified 50th-percentile bar chart, restricted to this stratum.",
                how_to_read="Same layout as the arm's unstratified bar chart.",
            ))
            blocks.append(bar50)
            blocks.append(plot_result(
                f"Mean knockdown for {primary_label} vs. {control_label} within the {stratum_label} stratum "
                f"is {_kd50:+.1f} percentage points (n={_strat_n:,} {primary_label} cells)."
            ))
        blocks.append(_tidy_strat_table(_strat_sub))
        return mo.vstack(blocks)

    _STRATUM_ORDER = ["All cells (unstratified)", "ASGR1+", "ALB+", "ASGR1+/ALB+"]

    def _arm_strat_comparison(cfg, plot_n):
        """Single at-a-glance grouped bar chart: guide-active group's knockdown
        (1st pctile red, 50th pctile purple) across every stratum that applies
        to this arm."""
        primary_label = cfg["groups"][-1][0]
        _sub = marker_strat_df[(marker_strat_df["arm"] == cfg["arm_label"]) & (marker_strat_df["label"] == primary_label)]
        agg = _sub.groupby("stratum", as_index=False)[["knockdown_pp_1st", "knockdown_pp_50th"]].mean()
        agg["stratum"] = pd.Categorical(agg["stratum"], categories=_STRATUM_ORDER, ordered=True)
        agg = agg.sort_values("stratum")

        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=agg["stratum"].astype(str), y=agg["knockdown_pp_1st"], name="1st pctile (tail effect)",
            marker_color="red", text=[f"{v:.1f}pp" for v in agg["knockdown_pp_1st"]], textposition="outside",
        ))
        fig.add_trace(go.Bar(
            x=agg["stratum"].astype(str), y=agg["knockdown_pp_50th"], name="50th pctile (population shift)",
            marker_color="purple", text=[f"{v:.1f}pp" for v in agg["knockdown_pp_50th"]], textposition="outside",
        ))
        _all_vals = list(agg["knockdown_pp_1st"]) + list(agg["knockdown_pp_50th"])
        _all_vals = [v for v in _all_vals if v == v]
        ymax = max(50, max(_all_vals) * 1.25) if _all_vals else 50
        ymin = min(0, min(_all_vals) * 1.25) if _all_vals else 0
        fig.update_layout(
            title=f"Plot {plot_n}. {cfg['arm_label']}: knockdown across strata ({primary_label} vs. mock)",
            barmode="group", yaxis_title="Knockdown (pp)", yaxis_range=[ymin, ymax],
            height=380, margin=dict(t=60),
            legend=dict(itemclick="toggle", itemdoubleclick="toggleothers"),
        )
        return fig, agg

    def _arm_strat_tldr(cfg, agg):
        _base = agg[agg["stratum"] == "All cells (unstratified)"]
        if not len(_base):
            return mo.md("**TL;DR:** insufficient data to compare strata.")
        base1 = float(_base["knockdown_pp_1st"].iloc[0])
        base50 = float(_base["knockdown_pp_50th"].iloc[0])
        parts = []
        for strat in ["ASGR1+", "ALB+", "ASGR1+/ALB+"]:
            row = agg[agg["stratum"] == strat]
            if not len(row):
                continue
            v1 = float(row["knockdown_pp_1st"].iloc[0])
            v50 = float(row["knockdown_pp_50th"].iloc[0])
            d1, d50 = v1 - base1, v50 - base50
            dir1 = "lowers" if d1 < -0.5 else "raises" if d1 > 0.5 else "leaves essentially unchanged"
            dir50 = "lowers" if d50 < -0.5 else "raises" if d50 > 0.5 else "leaves essentially unchanged"
            parts.append(
                f"**{strat}** selection {dir1} 1st-percentile knockdown ({base1:.1f}pp -> {v1:.1f}pp, {d1:+.1f}pp) "
                f"and {dir50} the 50th-percentile metric ({base50:.1f}pp -> {v50:.1f}pp, {d50:+.1f}pp)."
            )
        return mo.md("**TL;DR:** " + " ".join(parts))

    def build_arm_with_strat(cfg, base_content):
        """Wrap an arm's existing unstratified content ('All cells') together
        with ASGR1+/ALB+/double-positive tabs, where valid data exists -- folded
        into the arm's own tab rather than a separate standalone section."""
        asgr1_gate, alb_gate = cfg["asgr1_gate"], cfg["alb_gate"]
        tabs = {"All cells": base_content}
        running = cfg["plot_start"]
        if asgr1_gate is not None:
            tabs["ASGR1+"] = build_marker_stratum_block(cfg, "ASGR1+", ("PE-A", asgr1_gate), None, running)
            running += 4
        if alb_gate is not None:
            tabs["ALB+"] = build_marker_stratum_block(cfg, "ALB+", None, ("APC-A", alb_gate), running)
            running += 4
        if asgr1_gate is not None and alb_gate is not None:
            tabs["ASGR1+/ALB+"] = build_marker_stratum_block(cfg, "ASGR1+/ALB+", ("PE-A", asgr1_gate), ("APC-A", alb_gate), running)
            running += 4
        if len(tabs) == 1:
            return base_content

        comparison_fig, agg = _arm_strat_comparison(cfg, running)
        tldr = _arm_strat_tldr(cfg, agg)
        _primary_label = cfg["groups"][-1][0]

        return mo.vstack([
            plot_overview(
                running, f"{cfg['arm_label']}: knockdown across strata",
                why=f"States whether selecting for a hepatocyte-marker-positive cell group changes the measured knockdown for {cfg['arm_label']}, compared to using all gated cells.",
                how=f"The 1st- and 50th-percentile-of-mock knockdown (pp) for the {_primary_label} group is computed once per stratum (all cells, and each marker-positive subset), then plotted side by side.",
                how_to_read="Each group of bars is one cell stratum; red bars are the 1st-percentile-of-mock metric, purple bars are the 50th-percentile-of-mock metric, both in percentage points of knockdown.",
            ),
            comparison_fig,
            tldr,
            mo.md(
                "*ASGR1-FITC is excluded (non-specific clone) -- only the alternate "
                "ASGR1-PE channel and/or Albumin-AF647 are used for the "
                "marker-stratified tabs below, each AND'd with this arm's own "
                "infection/guide-delivery gate.*"
            ),
            mo.ui.tabs(tabs),
        ])


    return (build_arm_with_strat,)


@app.cell(hide_code=True)
def _(
    STRAT_ARM_CFGS,
    arm1_content_base,
    arm2_content_base,
    arm3_content_base,
    arm4_content_base,
    arm6_content_base,
    build_arm_with_strat,
    mo,
):
    _cfg_by_label = {c["arm_label"]: c for c in STRAT_ARM_CFGS}
    arm1_content = build_arm_with_strat(_cfg_by_label["Arm 1: pDRT103/106, 8_3"], arm1_content_base)
    arm2_content = build_arm_with_strat(_cfg_by_label["Arm 2: AA173+AA239, 17_3"], arm2_content_base)
    arm3_content = build_arm_with_strat(_cfg_by_label["Arm 3: AA173+AA239, 8_3"], arm3_content_base)
    # Arm 4 (G4-G7) has real, verified Albumin-AF647 staining -- stratified by
    # ALB+ (no valid ASGR1 channel available, FITC-only clone excluded).
    # Arm 5 stays unstratified: no valid marker data at all on H3-H8.
    arm4_content = build_arm_with_strat(_cfg_by_label["Arm 4: AA228 mCherry, 8_3"], arm4_content_base)
    arm6_content = build_arm_with_strat(_cfg_by_label["Arm 6: WTC11 fully-integrated, C1-C4"], arm6_content_base)
    mo.md("Marker-stratified content folded into arms 1, 2, 3, 6.")
    return arm1_content, arm2_content, arm3_content, arm4_content, arm6_content


@app.cell(hide_code=True)
def _(debris_slider, go, mo, np, plot_overview, plot_result, raw_wells):
    _pooled_fsc = np.concatenate([d["FSC-A"].values for d in raw_wells.values()])
    _pooled_fsc = _pooled_fsc[(_pooled_fsc > -50_000) & (_pooled_fsc < 3_000_000)]
    _fig = go.Figure()
    _fig.add_trace(go.Histogram(x=_pooled_fsc, nbinsx=150, marker_color="#4C78A8"))
    _fig.add_shape(
        type="line", x0=debris_slider.value, x1=debris_slider.value, y0=0, y1=1,
        yref="paper", line=dict(color="red", width=3, dash="dash"), layer="above",
    )
    _fig.add_annotation(
        x=debris_slider.value, y=1.03, yref="paper",
        text=f"debris gate = {debris_slider.value:,.0f}", showarrow=False,
        font=dict(color="red"),
    )
    _fig.update_layout(
        title="Plot 14. Pooled FSC-A -- debris/cell valley check",
        xaxis_title="FSC-A (linear)", yaxis_title="count", height=380, margin=dict(t=60),
    )
    _pct_below_debris = 100 * float(np.mean(_pooled_fsc <= debris_slider.value))
    debris_check = mo.vstack([
        plot_overview(
            14, "Pooled FSC-A -- debris/cell valley check",
            why="Confirms the FSC-A debris/cell floor sits in an actual valley of the pooled distribution before trusting any downstream gate.",
            how="FSC-A values from every loaded well are pooled and histogrammed; the dashed red line marks the current debris_slider threshold.",
            how_to_read="The valley on FSC-A sits roughly 200,000-300,000 on this dataset's linear scale; the gate should sit in that valley, not on either peak.",
        ),
        _fig,
        plot_result(f"{_pct_below_debris:.1f}% of pooled events fall at or below the current debris gate (excluded as debris)."),
    ])
    debris_check
    return


@app.cell(hide_code=True)
def _(add_threshold, debris_slider, go, mo, np, plot_overview, raw_wells):
    _rng_d = np.random.default_rng(0)
    _xs_d, _ys_d = [], []
    for _d in raw_wells.values():
        _n = len(_d)
        _cap = min(_n, 2500)
        _idx = _rng_d.choice(_n, _cap, replace=False)
        _xs_d.append(_d["FSC-A"].values[_idx])
        _ys_d.append(_d["SSC-A"].values[_idx])
    _xs_d = np.concatenate(_xs_d)
    _ys_d = np.concatenate(_ys_d)
    _fig_ds = go.Figure()
    _fig_ds.add_trace(go.Scattergl(
        x=_xs_d, y=_ys_d, mode="markers",
        marker=dict(size=2, opacity=0.25, color="#4C78A8"),
        name="events",
    ))
    add_threshold(_fig_ds, debris_slider.value, f"FSC-A debris gate = {debris_slider.value:,.0f}")
    _fig_ds.add_shape(
        type="line", x0=0, x1=1, xref="paper", y0=0, y1=0,
        line=dict(color="orange", width=2, dash="dot"), layer="above",
    )
    _fig_ds.add_annotation(
        x=0.02, y=0, xref="paper", text="SSC-A = 0 floor", showarrow=False,
        font=dict(color="orange", size=10), xanchor="left", yanchor="bottom",
    )
    _fig_ds.update_layout(
        title="Plot 15. Pooled FSC-A vs SSC-A -- debris/cell valley in 2D",
        xaxis_title="FSC-A", yaxis_title="SSC-A", height=420, margin=dict(t=60),
    )
    _fig_ds.add_annotation(
        text="downsampled to <=2,500 events/well for display", x=0, y=-0.14,
        xref="paper", yref="paper", showarrow=False, font=dict(size=10, color="gray"),
    )
    debris_scatter = mo.vstack([
        plot_overview(
            15, "Pooled FSC-A vs SSC-A -- debris/cell valley in 2D",
            why="Checks the same debris/cell separation as Plot 14, but in 2D against SSC-A, to catch any population the 1D FSC-A view alone would miss.",
            how="A downsampled pool of events from every loaded well, plotted as FSC-A vs SSC-A, with the debris gate and SSC-A=0 floor marked.",
            how_to_read="The dashed red vertical line is the FSC-A debris gate; the dotted orange horizontal line is the SSC-A=0 floor.",
        ),
        _fig_ds,
    ])
    debris_scatter
    return


@app.cell(hide_code=True)
def _(mo):
    debris_slider = mo.ui.slider(
        start=0, stop=1_500_000, step=10_000, value=250_000,
        label="FSC-A debris/cell gate (events at/below this FSC-A value are excluded as debris)",
        full_width=True, show_value=True,
    )
    debris_slider
    return (debris_slider,)


@app.cell(hide_code=True)
def _(mo):
    doublet_slider = mo.ui.slider(
        start=0, stop=8000, step=100, value=3800,
        label="FSC-Width doublet/clump cutoff (events above this FSC-Width are excluded as doublets)",
        full_width=True, show_value=True,
    )
    doublet_slider
    return (doublet_slider,)


@app.cell(hide_code=True)
def _(doublet_slider, go, mo, np, plot_overview, plot_result, raw_wells):
    _pooled_w = np.concatenate([d["FSC-Width"].values for d in raw_wells.values()])
    _pooled_w = _pooled_w[(_pooled_w > 0) & (_pooled_w < 8000)]
    _fig_w = go.Figure()
    _fig_w.add_trace(go.Histogram(x=_pooled_w, nbinsx=150, marker_color="#54A24B"))
    _fig_w.add_shape(
        type="line", x0=doublet_slider.value, x1=doublet_slider.value, y0=0, y1=1,
        yref="paper", line=dict(color="red", width=3, dash="dash"), layer="above",
    )
    _fig_w.add_annotation(
        x=doublet_slider.value, y=1.03, yref="paper",
        text=f"doublet cutoff = {doublet_slider.value:,.0f}", showarrow=False,
        font=dict(color="red"),
    )
    _fig_w.update_layout(
        title="Plot 16. Pooled FSC-Width -- doublet/clump cutoff check",
        xaxis_title="FSC-Width (linear)", yaxis_title="count", height=380, margin=dict(t=60),
    )
    _pct_above_doublet = 100 * float(np.mean(_pooled_w > doublet_slider.value))
    doublet_check = mo.vstack([
        plot_overview(
            16, "Pooled FSC-Width -- doublet/clump cutoff check",
            why="Confirms the FSC-Width doublet/clump cutoff clips the distribution's tail rather than cutting into the main singlet peak.",
            how="FSC-Width values from every loaded well are pooled and histogrammed; the dashed red line marks the current doublet_slider threshold.",
            how_to_read="The pooled distribution has one main peak (~900-1000, true singlets) and a smoothly decaying right tail; the cutoff should sit past the main peak.",
        ),
        _fig_w,
        plot_result(f"{_pct_above_doublet:.1f}% of pooled events exceed the current FSC-Width cutoff (excluded as doublets/clumps)."),
    ])
    doublet_check
    return


@app.cell(hide_code=True)
def _(doublet_slider, go, mo, np, plot_overview, raw_wells):
    _rng_w = np.random.default_rng(0)
    _xs_w, _ys_w = [], []
    for _d in raw_wells.values():
        _n = len(_d)
        _cap = min(_n, 2500)
        _idx = _rng_w.choice(_n, _cap, replace=False)
        _xs_w.append(_d["FSC-A"].values[_idx])
        _ys_w.append(_d["FSC-Width"].values[_idx])
    _xs_w = np.concatenate(_xs_w)
    _ys_w = np.concatenate(_ys_w)
    _fig_ws = go.Figure()
    _fig_ws.add_trace(go.Scattergl(
        x=_xs_w, y=_ys_w, mode="markers",
        marker=dict(size=2, opacity=0.25, color="#54A24B"),
        name="events",
    ))
    _fig_ws.add_shape(
        type="line", x0=0, x1=1, xref="paper", y0=doublet_slider.value, y1=doublet_slider.value,
        line=dict(color="red", width=3, dash="dash"), layer="above",
    )
    _fig_ws.add_annotation(
        x=0.02, y=doublet_slider.value, xref="paper",
        text=f"FSC-Width doublet cutoff = {doublet_slider.value:,.0f}", showarrow=False,
        font=dict(color="red"), xanchor="left", yanchor="bottom",
    )
    _fig_ws.update_layout(
        title="Plot 17. Pooled FSC-A vs FSC-Width -- doublet band above the cutoff",
        xaxis_title="FSC-A", yaxis_title="FSC-Width", height=420, margin=dict(t=60),
    )
    _fig_ws.add_annotation(
        text="downsampled to <=2,500 events/well for display", x=0, y=-0.14,
        xref="paper", yref="paper", showarrow=False, font=dict(size=10, color="gray"),
    )
    doublet_scatter = mo.vstack([
        plot_overview(
            17, "Pooled FSC-A vs FSC-Width -- doublet band above the cutoff",
            why="Checks the same doublet cutoff as Plot 16, but in 2D against FSC-A, to confirm the doublet band sits consistently above the cutoff across the FSC-A range.",
            how="A downsampled pool of events from every loaded well, plotted as FSC-A vs FSC-Width, with the doublet cutoff marked.",
            how_to_read="The dashed red horizontal line is the FSC-Width doublet cutoff; events above it are excluded as doublets/clumps.",
        ),
        _fig_ws,
    ])
    doublet_scatter
    return


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import numpy as np
    import pandas as pd
    import requests
    import urllib.parse
    import io
    import plotly.graph_objects as go

    return go, io, mo, np, pd, requests, urllib


@app.cell(hide_code=True)
def _(mo):
    token_input = mo.ui.text(
        kind="password",
        label="GCS access token (paste output of `gcloud auth print-access-token`; expires ~1hr, re-paste if you hit 401s)",
        full_width=True,
    )
    token_input
    return (token_input,)


@app.cell(hide_code=True)
def _(requests, token_input, urllib):
    BUCKET = "20240617-landerlab-crispri-project"
    PREFIX = "Flow_Data_CRISPRi/20260928_HLC_Trial3_D27_fixed/"

    def gcs_get(name: str) -> bytes:
        enc = urllib.parse.quote(name, safe="")
        url = f"https://storage.googleapis.com/storage/v1/b/{BUCKET}/o/{enc}?alt=media"
        headers = {"Authorization": f"Bearer {token_input.value}"}
        r = requests.get(url, headers=headers)
        r.raise_for_status()
        return r.content

    return PREFIX, gcs_get


@app.cell(hide_code=True)
def _(np, pd):
    def parse_fcs_bytes(data: bytes):
        text_start = int(data[10:18].decode().strip())
        text_end = int(data[18:26].decode().strip())
        data_start_hdr = int(data[26:34].decode().strip())
        data_end_hdr = int(data[34:42].decode().strip())
        text = data[text_start:text_end + 1].decode("latin-1")
        delim = text[0]
        parts = text[1:].split(delim)
        kv = {}
        it = iter(parts)
        for k, v in zip(it, it):
            kv[k.upper()] = v
        data_start = data_start_hdr or int(kv["$BEGINDATA"])
        data_end = data_end_hdr or int(kv["$ENDDATA"])
        npar = int(kv["$PAR"]); ntot = int(kv["$TOT"])
        endian = "<" if kv["$BYTEORD"].startswith("1,2") else ">"
        dt_map = {"F": "f4", "D": "f8", "I": "u4"}
        dtype = np.dtype(endian + dt_map[kv["$DATATYPE"]])
        raw = np.frombuffer(data[data_start:data_end + 1], dtype=dtype).reshape(ntot, npar)
        names = [kv.get(f"$P{i}S") or kv.get(f"$P{i}N", f"P{i}") for i in range(1, npar + 1)]
        return pd.DataFrame(raw, columns=names), kv

    return (parse_fcs_bytes,)


@app.cell(hide_code=True)
def _(PREFIX, gcs_get, io, pd):
    _xls_bytes = gcs_get(PREFIX + "20260928_fixed_SampleID.xls")
    sample_sheet = pd.read_excel(io.BytesIO(_xls_bytes), engine="xlrd", header=0)
    sample_sheet["Well"] = sample_sheet["Well ID"].str.replace("03-Well-", "", regex=False)
    sample_sheet = sample_sheet.set_index("Well", drop=False)
    sample_sheet
    return (sample_sheet,)


@app.cell(hide_code=True)
def _(PREFIX, gcs_get, parse_fcs_bytes):
    # Only the wells needed for the CD81 knockdown arms + fully-unstained naive
    # references per cell line (A1=WTC11, A2=17_3, A3=8_3) + D1-D4 (single-stain
    # compensation controls, 8_3: D1=CD81-BV421, D2=ASGR1-FITC, D3=ASGR1-PE,
    # D4=Albumin-AF647/read-in-APC-A) + D5 (Thy1.1-FITC stain, no virus, 8_3) as
    # an extra background-check well + A4/A8/D5 (Thy1.1-FITC stained, no virus --
    # FMO-style stained-but-uninfected references for the Thy1.1 gate, per cell
    # line: A4=WTC11, A8=17_3, D5=8_3).
    # Explicitly excludes: the unfixed sibling folder, K562 wells (F1-F4), and all
    # other ASGR1/Albumin-only wells not needed for compensation -- out of scope
    # for this CD81-only analysis.
    WELLS = [
        "A1", "A2", "A3", "A4", "A8",
        "B1", "B2", "B3", "B4", "B5", "B6",
        "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8",
        "C9", "C10", "C11", "C12",
        "D1", "D2", "D3", "D4", "D5",
        "G4", "G5", "G6", "G7",
        "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8",
    ]

    def _load_well(w):
        _data = gcs_get(PREFIX + f"03-Well-{w}.fcs")
        _df, _kv = parse_fcs_bytes(_data)
        return _df

    raw_wells = {w: _load_well(w) for w in WELLS}
    f"Loaded {len(raw_wells)} wells, {sum(len(d) for d in raw_wells.values()):,} total events"
    return WELLS, raw_wells


@app.cell(hide_code=True)
def _(
    WELLS,
    debris_slider,
    doublet_slider,
    go,
    mo,
    np,
    pd,
    raw_wells,
    sample_sheet,
    ssc_cap_slider,
):
    MIN_CELLS_PER_GATE = 300
    _CONTROL_TOKENS = ["ork", "non-targeting", "nt-", "unstained", "no guide", "parental", "control"]

    def is_control_label(text: str) -> bool:
        t = (text or "").lower()
        return any(tok in t for tok in _CONTROL_TOKENS)

    # --- Compensation: least-squares spillover matrix from single-stain wells ---
    # D1-D4 (8_3, fixation-matched) are each stained with exactly one fluorophore
    # used as a CD81/ASGR1/Albumin readout channel elsewhere in this notebook.
    # Albumin's Alexa Fluor 647 secondary reads almost entirely in the APC-A
    # channel on this instrument (confirmed empirically: well A7/D4 shows no
    # bimodal signal anywhere except APC-A), so D4 doubles as the APC-A
    # single-stain reference -- no separate APC-only control well exists on this
    # plate, but this one is equivalent for spillover-estimation purposes.
    _COMP_CHANNELS = ["BV421-A", "FITC-A", "PE-A", "APC-A"]
    _SINGLE_STAIN_WELLS = {"BV421-A": "D1", "FITC-A": "D2", "PE-A": "D3", "APC-A": "D4"}
    _COMP_NAIVE_WELL = "A3"

    def _build_spillover_matrix():
        naive = raw_wells[_COMP_NAIVE_WELL]
        naive_medians = {ch: float(np.median(naive[ch].values)) for ch in _COMP_CHANNELS}
        n = len(_COMP_CHANNELS)
        S = np.eye(n)
        for i, target_ch in enumerate(_COMP_CHANNELS):
            stained = raw_wells[_SINGLE_STAIN_WELLS[target_ch]]
            primary_signal = stained[target_ch].values - naive_medians[target_ch]
            primary_signal = np.clip(primary_signal, 1.0, None)
            for j, other_ch in enumerate(_COMP_CHANNELS):
                if i == j:
                    continue
                other_signal = stained[other_ch].values - naive_medians[other_ch]
                slope = float(np.polyfit(primary_signal, other_signal, 1)[0])
                S[j, i] = max(0.0, slope)
        return S

    SPILLOVER_MATRIX = pd.DataFrame(0.0, index=_COMP_CHANNELS, columns=_COMP_CHANNELS)
    _SPILLOVER_RAW = _build_spillover_matrix()
    for _i, _ch in enumerate(_COMP_CHANNELS):
        SPILLOVER_MATRIX.iloc[_i] = _SPILLOVER_RAW[_i]
    _INV_SPILLOVER = np.linalg.inv(_SPILLOVER_RAW)

    compensation_table = mo.ui.table(SPILLOVER_MATRIX.round(3), selection=None)

    def compensate(d: pd.DataFrame) -> pd.DataFrame:
        """Apply the inverse spillover matrix to the 4 compensated channels;
        every other channel (FSC/SSC/Time/other fluorophores not in this panel's
        4-channel readout set) passes through unchanged."""
        d = d.copy()
        raw = d[_COMP_CHANNELS].values
        comp = raw @ _INV_SPILLOVER.T
        for i, ch in enumerate(_COMP_CHANNELS):
            d[ch] = comp[:, i]
        return d

    # --- Singlet gate: FSC-H:FSC-A ratio band, calibrated empirically ---
    def _pre_singlet_gate(well: str) -> pd.DataFrame:
        d = compensate(raw_wells[well])
        mask = (
            (d["FSC-A"].values > debris_slider.value)
            & (d["SSC-A"].values > 0)
            & (d["SSC-A"].values <= ssc_cap_slider.value)
            & (d["FSC-Width"].values <= doublet_slider.value)
        )
        return d[mask]

    def _singlet_ratio_band():
        ratios = []
        for w in WELLS:
            d = _pre_singlet_gate(w)
            if len(d):
                ratios.append(d["FSC-H"].values / np.maximum(d["FSC-A"].values, 1.0))
        pooled = np.concatenate(ratios) if ratios else np.array([1.0])
        return float(np.percentile(pooled, 2.5)), float(np.percentile(pooled, 97.5))

    SINGLET_RATIO_BAND = _singlet_ratio_band()

    def debris_gate(well: str) -> pd.DataFrame:
        """Cell gate (FSC-A floor + SSC-A upper cap, a rectangle not just a
        floor) + doublet gate (FSC-Width ceiling) + singlet gate (FSC-H:FSC-A
        ratio band, central ~95% of the pooled ratio across all loaded wells),
        applied to compensated channel values."""
        d = _pre_singlet_gate(well)
        if len(d):
            ratio = d["FSC-H"].values / np.maximum(d["FSC-A"].values, 1.0)
            d = d[(ratio >= SINGLET_RATIO_BAND[0]) & (ratio <= SINGLET_RATIO_BAND[1])]
        return d

    def pctile_gate(well: str, chan: str, pct: float = 99.0) -> float:
        d = debris_gate(well)
        return float(np.percentile(d[chan].values, pct))

    def apply_markers(well: str, marker_specs: list[tuple[str, float]]) -> pd.DataFrame:
        d = debris_gate(well)
        mask = np.ones(len(d), dtype=bool)
        for chan, thr in marker_specs:
            mask &= d[chan].values > thr
        return d[mask]

    def gated_range(naive_vals, pooled_vals, pad_frac: float = 0.10, hi_pct: float = 99.0, lo_pct: float = 1.0):
        """Right edge: hi_pct-of-pooled-data plus padding (existing convention).
        Left edge: lo_pct (default 1st percentile) of the UNSTAINED control's own
        signal, no further left padding -- keeps the plot from being dominated by
        the unstained reference's low-signal tail."""
        naive_vals = np.asarray(naive_vals)
        pooled_vals = np.asarray(pooled_vals)
        lo_anchor = np.percentile(naive_vals, lo_pct) if len(naive_vals) else np.percentile(pooled_vals, lo_pct)
        hi_anchor = np.percentile(pooled_vals, hi_pct)
        span = max(hi_anchor - lo_anchor, 1.0)
        return lo_anchor, hi_anchor + pad_frac * span

    def add_threshold(fig, x, label, color="red", y=1.05):
        fig.add_shape(
            type="line", x0=x, x1=x, y0=0, y1=1, yref="paper",
            line=dict(color=color, width=3, dash="dash"), layer="above",
        )
        fig.add_annotation(
            x=x, y=y, yref="paper", text=label, showarrow=False, font=dict(color=color),
        )

    def interactive_hist(
        traces: dict, gate: float | None, gate_label: str, xrange, title: str, nbins: int = 60,
        gate2: float | None = None, gate2_label: str = "", gate2_color: str = "purple",
    ):
        fig = go.Figure()
        colors = ["#4C78A8", "#E45756", "#54A24B", "#F58518", "#B279A2", "#72B7B2"]
        for i, (label, vals) in enumerate(traces.items()):
            vals = np.asarray(vals)
            vals = vals[(vals >= xrange[0]) & (vals <= xrange[1])]
            counts, edges = np.histogram(vals, bins=nbins, range=xrange)
            counts = counts / counts.max() * 100 if counts.max() > 0 else counts
            centers = (edges[:-1] + edges[1:]) / 2
            fig.add_trace(go.Scatter(
                x=centers, y=counts, mode="lines", fill="tozeroy",
                name=f"{label} (n={len(vals):,})", line=dict(color=colors[i % len(colors)]),
                opacity=0.55,
            ))
        if gate is not None:
            add_threshold(fig, gate, gate_label, color="red", y=1.05)
        if gate2 is not None:
            add_threshold(fig, gate2, gate2_label, color=gate2_color, y=1.12)
        fig.update_layout(
            title=title, xaxis_title="signal", yaxis_title="% of peak (normalized)",
            height=380, margin=dict(t=70), xaxis_range=list(xrange),
            legend=dict(itemclick="toggle", itemdoubleclick="toggleothers"),
        )
        return fig

    def scatter_gate(
        traces: dict, xrange, title: str, y_chan: str = "SSC-A",
        gate: float | None = None, gate_label: str = "",
        gate2: float | None = None, gate2_label: str = "", gate2_color: str = "purple",
    ):
        """Density-contour scatter: one Histogram2dContour per pooled group
        (lines-only, distinct colors) instead of overplotted raw markers.
        Y-axis capped at the 99.5th percentile of the plotted events (display
        only) so rare high-SSC-A outliers don't compress the main population."""
        fig = go.Figure()
        colors = ["#4C78A8", "#E45756", "#54A24B", "#F58518", "#B279A2", "#72B7B2"]
        all_y = np.concatenate([np.asarray(yv) for _, yv in traces.values()]) if traces else np.array([0.0])
        y_cap = float(np.percentile(all_y, 99.5)) if len(all_y) else 1.0
        for i, (label, (xv, yv)) in enumerate(traces.items()):
            xv = np.asarray(xv); yv = np.asarray(yv)
            if not len(xv):
                continue
            _color = colors[i % len(colors)]
            fig.add_trace(go.Histogram2dContour(
                x=xv, y=yv, name=f"{label} (n={len(xv):,})",
                contours=dict(coloring="lines", showlines=True),
                line=dict(width=2, color=_color), colorscale=[[0, _color], [1, _color]],
                showscale=False, showlegend=True, ncontours=8,
            ))
        if gate is not None:
            add_threshold(fig, gate, gate_label, color="red", y=1.05)
        if gate2 is not None:
            add_threshold(fig, gate2, gate2_label, color=gate2_color, y=1.12)
        fig.update_layout(
            title=title, xaxis_title="signal", yaxis_title=y_chan, height=420,
            margin=dict(t=70), xaxis_range=list(xrange), yaxis_range=[0, y_cap],
            legend=dict(itemclick="toggle", itemdoubleclick="toggleothers"),
        )
        return fig

    def tidy_summary_table(df: pd.DataFrame):
        _df = df.copy()
        for c in _df.columns:
            if "pct" in c.lower() or c == "mfi" or c == "mfi_ratio_vs_control":
                _df[c] = _df[c].round(1)
        rename = {
            "well": "Well", "label": "Condition", "n": "n (gated cells)",
            "pct_below_gate": "% below 1st-pctile-mock gate",
            "knockdown_pp": "Knockdown, 1st pctile (pp vs mock)",
            "pct_below_gate50": "% below 50th-pctile-mock gate",
            "knockdown_pp_50": "Knockdown, 50th pctile (pp vs mock)",
            "mfi": "Median readout (MFI)",
            "mfi_ratio_vs_control": "MFI ratio vs control",
        }
        _df = _df.rename(columns={k: v for k, v in rename.items() if k in _df.columns})
        return mo.ui.table(_df, selection=None)

    def replicate_labels(wells: list[str], group_label: str = None) -> dict:
        """Reader-friendly per-replicate labels. When group_label is given (the
        normal case from build_arm), labels never include a well ID -- just
        "{group_label} rep{i}". Falls back to the legacy cell-line+well format
        only when no group label is supplied."""
        labels = {}
        for i, w in enumerate(wells, start=1):
            if group_label:
                labels[w] = f"{group_label} rep{i}"
            else:
                cell_line = sample_sheet.loc[w, "Cell Line"]
                labels[w] = f"{cell_line} rep{i} ({w})"
        return labels

    _CHANNEL_COFACTORS = {}

    def channel_cofactor(channel: str) -> float:
        """One fixed biexponential cofactor per channel, derived from that
        channel's signal in the shared unstained reference well (A3), cached and
        reused across every plot of that channel -- not recomputed per-plot from
        that plot's own gate value."""
        if channel not in _CHANNEL_COFACTORS:
            naive_vals = debris_gate(_COMP_NAIVE_WELL)[channel].values
            ref = max(float(np.percentile(naive_vals, 95)), 1.0) if len(naive_vals) else 1.0
            _CHANNEL_COFACTORS[channel] = ref / 5.0
        return _CHANNEL_COFACTORS[channel]

    def biexp(x, cofactor: float):
        """arcsinh-based biexponential display transform: linear-like near zero,
        log-like at high magnitude, handles negative/near-zero compensated
        values gracefully (unlike a true log axis)."""
        return np.arcsinh(np.asarray(x, dtype=float) / cofactor)

    def biexp_ticks(lo_raw: float, hi_raw: float, cofactor: float):
        """Geometric raw-value tick marks (0, cofactor, 3*cofactor, 9*cofactor, ...)
        re-expressed at their transformed positions, so axis labels show real
        fluorescence units instead of meaningless arcsinh-space numbers."""
        vals = [0.0]
        mult = 1.0
        while True:
            v = cofactor * mult
            if v > hi_raw * 1.02:
                break
            vals.append(v)
            mult *= 3.0
        if lo_raw < -cofactor * 0.5:
            vals = [lo_raw] + vals
        vals = sorted(set(round(v, 0) for v in vals))
        tickvals = [biexp(v, cofactor) for v in vals]
        ticktext = [f"{v:,.0f}" for v in vals]
        return tickvals, ticktext

    def _condition_label(well: str) -> str:
        """Human-readable construct/condition label for a well -- cell line +
        virus/construct + integrated-guide + Dox status -- reused anywhere a bare
        well ID would otherwise show up in a chart legend."""
        row = sample_sheet.loc[well]
        cell_line = str(row["Cell Line"])

        def _clean(x):
            s = str(x).strip() if x is not None else ""
            return s if s and s.lower() not in ("nan", "no", "none") else None

        parts = [cell_line]
        v1c, v2c = _clean(row.get("Virus 1 Type")), _clean(row.get("Virus 2 Type"))
        if v1c:
            parts.append(v1c.replace(" ", "-"))
        if v2c:
            parts.append(v2c.replace(" ", "-"))
        guide_c = _clean(row.get("Integrated Guide"))
        if guide_c:
            parts.append(f"integrated-{guide_c}")
        label = " + ".join(parts)
        if _clean(row.get("Doxycycline Induction")) == "Yes":
            label += ", Dox-on"
        return label

    def calibration_trace_labels(wells: list[str]) -> dict:
        """Map each well to a readable trace label (never a bare well ID),
        appending '(repN)' only when multiple wells in this list share the
        identical condition label."""
        base = {w: _condition_label(w) for w in wells}
        totals = {}
        for w in wells:
            totals[base[w]] = totals.get(base[w], 0) + 1
        seen = {}
        out = {}
        for w in wells:
            b = base[w]
            if totals[b] > 1:
                seen[b] = seen.get(b, 0) + 1
                out[w] = f"{b} (rep{seen[b]})"
            else:
                out[w] = b
        return out

    def plot_overview(n, title, why, how, how_to_read):
        return mo.md(
            f"**Plot {n}. {title}**\n\n"
            f"*Why:* {why}\n\n"
            f"*How:* {how}\n\n"
            f"*How to read:* {how_to_read}"
        )

    def plot_result(text):
        return mo.md(f"*Result:* {text}")

    def gating_hierarchy_table(wells: list[str]):
        """Per well: event counts and % of the previous step through all events
        -> cells (FSC-A floor + SSC-A cap) -> single cells (+ doublet width cut
        + FSC-H:FSC-A singlet band)."""
        rows = []
        for w in wells:
            n_all = len(raw_wells[w])
            n_cell = len(_pre_singlet_gate(w))
            n_single = len(debris_gate(w))
            rows.append({
                "Well": w,
                "All events": n_all,
                "Cells (FSC-A floor + SSC-A cap)": n_cell,
                "% of all events": round(100 * n_cell / n_all, 1) if n_all else 0.0,
                "Single cells (+ width cut + FSC-H:FSC-A band)": n_single,
                "% of cells": round(100 * n_single / n_cell, 1) if n_cell else 0.0,
            })
        return mo.ui.table(pd.DataFrame(rows), selection=None)


    return (
        MIN_CELLS_PER_GATE,
        add_threshold,
        apply_markers,
        biexp,
        biexp_ticks,
        calibration_trace_labels,
        channel_cofactor,
        debris_gate,
        gated_range,
        gating_hierarchy_table,
        interactive_hist,
        is_control_label,
        pctile_gate,
        plot_overview,
        plot_result,
        replicate_labels,
        scatter_gate,
        tidy_summary_table,
    )


@app.cell(hide_code=True)
def _(
    MIN_CELLS_PER_GATE,
    apply_markers,
    biexp,
    biexp_ticks,
    channel_cofactor,
    debris_gate,
    gated_range,
    gating_hierarchy_table,
    go,
    interactive_hist,
    is_control_label,
    mo,
    np,
    pd,
    plot_overview,
    plot_result,
    replicate_labels,
    sample_sheet,
    scatter_gate,
    tidy_summary_table,
):
    def build_arm(
        *,
        title: str,
        description: str,
        groups: list[tuple[str, list[str], list[tuple[str, float]]]],
        readout_channel: str,
        naive_well: str,
        control_group_index: int = 0,
        caveat: str | None = None,
        summary_flag: str = "",
        plot_start: int = 1,
    ):
        """groups: list of (group_label, wells, marker_specs) -- marker_specs is a
        list of (channel, threshold) applied as an AND gate (empty list = no
        infection gate, e.g. the fully-integrated stable line)."""
        naive_vals = debris_gate(naive_well)[readout_channel].values

        per_well_rows = []
        gated_by_well = {}
        for label, wells, marker_specs in groups:
            for w in wells:
                d = apply_markers(w, marker_specs) if marker_specs else debris_gate(w)
                gated_by_well[w] = (label, d)

        control_label, control_wells, control_marker_specs = groups[control_group_index]
        control_vals = np.concatenate([
            gated_by_well[w][1][readout_channel].values for w in control_wells
        ]) if all(len(gated_by_well[w][1]) > 0 for w in control_wells) else np.array([])
        gate = float(np.percentile(control_vals, 1)) if len(control_vals) else None
        gate50 = float(np.percentile(control_vals, 50)) if len(control_vals) else None
        control_mfi = float(np.median(control_vals)) if len(control_vals) else float("nan")

        pooled_vals = np.concatenate([d[readout_channel].values for _, d in gated_by_well.values() if len(d)])
        xr = gated_range(naive_vals, pooled_vals, hi_pct=99.0)

        # One fixed biexponential cofactor per channel (cached, derived from the
        # shared unstained reference well), not recomputed per-arm from this
        # arm's own gate value.
        cofactor = channel_cofactor(readout_channel)

        def _tx(v):
            return biexp(v, cofactor)

        naive_vals_t = _tx(naive_vals)
        xr_t = (float(_tx(xr[0])), float(_tx(xr[1])))
        gate_t = float(_tx(gate)) if gate is not None else None
        gate50_t = float(_tx(gate50)) if gate50 is not None else None

        # Histogram and scatter/density traces both pool replicate wells within
        # each condition group (one trace per group, not per well).
        hist_traces_by_group = {}
        scatter_traces_by_group = {}
        for w, (label, d) in gated_by_well.items():
            if len(d):
                vals_t = _tx(d[readout_channel].values)
                hist_traces_by_group.setdefault(label, []).append(vals_t)
                xs, ys = scatter_traces_by_group.setdefault(label, ([], []))
                xs.append(vals_t)
                ys.append(d["SSC-A"].values)
        hist_traces = {"unstained reference": naive_vals_t}
        hist_traces.update({label: np.concatenate(vs) for label, vs in hist_traces_by_group.items()})
        scatter_traces = {
            label: (np.concatenate(xs), np.concatenate(ys))
            for label, (xs, ys) in scatter_traces_by_group.items()
        }

        _gate1_label = f"1st %ile mock ({control_label}) = {gate:,.0f}" if gate else ""
        _gate50_label = f"50th %ile mock ({control_label}) = {gate50:,.0f}" if gate50 else ""

        fig_hist = interactive_hist(
            hist_traces, gate_t, _gate1_label,
            xr_t, f"Plot {plot_start}. {title}: {readout_channel} distribution",
            gate2=gate50_t, gate2_label=_gate50_label, gate2_color="purple",
        )
        fig_scatter = scatter_gate(
            scatter_traces, xr_t, f"Plot {plot_start + 1}. {title}: {readout_channel} vs SSC-A", y_chan="SSC-A",
            gate=gate_t, gate_label=_gate1_label, gate2=gate50_t, gate2_label=_gate50_label, gate2_color="purple",
        )

        _tickvals, _ticktext = biexp_ticks(xr[0], xr[1], cofactor)
        fig_hist.update_layout(xaxis=dict(
            tickvals=_tickvals, ticktext=_ticktext, title=f"{readout_channel} (biexponential scale, original units)",
        ))
        fig_scatter.update_layout(xaxis=dict(
            tickvals=_tickvals, ticktext=_ticktext, title=f"{readout_channel} (biexponential scale, original units)",
        ))

        excluded = []
        for label, wells, _ in groups:
            for w in wells:
                d = gated_by_well[w][1]
                n = len(d)
                if n < MIN_CELLS_PER_GATE and not is_control_label(label):
                    excluded.append((w, label, n))
                    continue
                vals = d[readout_channel].values
                pct = 100 * np.mean(vals < gate) if gate is not None and n else np.nan
                pct50 = 100 * np.mean(vals < gate50) if gate50 is not None and n else np.nan
                mfi = float(np.median(vals)) if n else float("nan")
                per_well_rows.append({
                    "well": w, "label": label, "n": n,
                    "pct_below_gate": pct, "pct_below_gate50": pct50,
                    "mfi": mfi, "mfi_ratio_vs_control": (mfi / control_mfi) if control_mfi and n else float("nan"),
                })

        df = pd.DataFrame(per_well_rows)
        if len(df) and gate is not None:
            ctrl_mean = df.loc[df["label"] == control_label, "pct_below_gate"].mean()
            df["knockdown_pp"] = df["pct_below_gate"] - ctrl_mean
        if len(df) and gate50 is not None:
            ctrl_mean50 = df.loc[df["label"] == control_label, "pct_below_gate50"].mean()
            df["knockdown_pp_50"] = df["pct_below_gate50"] - ctrl_mean50

        bar_fig = go.Figure()
        bar_fig_50 = go.Figure()
        if len(df):
            rl = {}
            for label, wells, _ in groups:
                rl.update(replicate_labels([w for w in wells if w in df["well"].values], group_label=label))
            x_labels = [rl.get(w, w) for w in df["well"]]
            colors_map = {g[0]: c for g, c in zip(groups, ["#4C78A8", "#E45756", "#54A24B", "#F58518"])}
            bar_fig.add_trace(go.Bar(
                x=x_labels, y=df["pct_below_gate"],
                text=[f"{v:.1f}%<br>(n={n:,})" for v, n in zip(df["pct_below_gate"], df["n"])],
                textposition="outside",
                marker_color=[colors_map.get(l, "#999") for l in df["label"]],
            ))
            ymax = max(50, float(df["pct_below_gate"].max()) * 1.25) if df["pct_below_gate"].notna().any() else 100
            bar_fig.update_layout(
                title=f"Plot {plot_start + 2}. {title}: % below 1st-percentile-of-mock gate (per replicate; tail effect)",
                yaxis_title="% below gate", yaxis_range=[0, min(100, ymax)],
                height=380, margin=dict(t=60),
            )
            bar_fig_50.add_trace(go.Bar(
                x=x_labels, y=df["pct_below_gate50"],
                text=[f"{v:.1f}%<br>(n={n:,})" for v, n in zip(df["pct_below_gate50"], df["n"])],
                textposition="outside",
                marker_color=[colors_map.get(l, "#999") for l in df["label"]],
            ))
            ymax50 = max(50, float(df["pct_below_gate50"].max()) * 1.25) if df["pct_below_gate50"].notna().any() else 100
            bar_fig_50.update_layout(
                title=f"Plot {plot_start + 3}. {title}: % below 50th-percentile-of-mock gate (per replicate; population-level shift)",
                yaxis_title="% below gate", yaxis_range=[0, min(100, ymax50)],
                height=380, margin=dict(t=60),
            )

        excl_note = ""
        if excluded:
            items = "; ".join(f"{lbl}, n={n}" for w, lbl, n in excluded)
            excl_note = f"*Excluded (population < {MIN_CELLS_PER_GATE} cells):* {items}"

        meta_cols = ["Well ID", "Sample ID #", "Cell Line", "Integrated Guide",
                     "Doxycycline Induction", "Virus 1 Type", "Virus 2 Type"]
        all_wells = [w for _, wells, _ in groups for w in wells]
        meta_table = mo.ui.table(sample_sheet.loc[all_wells, meta_cols], selection=None)
        hierarchy_table = gating_hierarchy_table(all_wells)

        _metric_note = mo.md(
            "*Two complementary metrics below: the red **1st-percentile-of-mock** "
            "gate captures a tail/low-signal effect (cells shifted furthest down); "
            "the purple **50th-percentile-of-mock** gate captures a broader "
            "population-level shift (whether the bulk of the distribution moved, "
            "not just the tail). Arms where the 50th-percentile knockdown is much "
            "larger than the 1st-percentile knockdown indicate a population-wide "
            "shift rather than a tail-only effect.*"
        )

        primary_label, primary_wells, _ = groups[-1]
        _kd1 = float(df.loc[df["label"] == primary_label, "knockdown_pp"].mean()) if len(df) and "knockdown_pp" in df else float("nan")
        _kd50 = float(df.loc[df["label"] == primary_label, "knockdown_pp_50"].mean()) if len(df) and "knockdown_pp_50" in df else float("nan")
        _n_primary = int(df.loc[df["label"] == primary_label, "n"].sum()) if len(df) else 0

        blocks = [mo.md(description)]
        if caveat:
            blocks.append(mo.callout(mo.md(caveat), kind="warn"))
        blocks.append(mo.ui.tabs({"Well metadata": meta_table, "Gating hierarchy": hierarchy_table}))
        blocks.append(_metric_note)
        blocks.append(plot_overview(
            plot_start, f"{title}: {readout_channel} distribution",
            why=f"States whether {primary_label} cells show a {readout_channel} shift relative to {control_label}, the arm's core knockdown question.",
            how=(
                f"{readout_channel} values (gated on this arm's guide/effector-delivery markers, singlet/cell-gated, and "
                f"compensated for cross-channel spillover) are pooled across replicate wells within each condition group, "
                f"then binned into a histogram and normalized so each group's tallest bin reads 100%. An unstained "
                f"reference well is overlaid for comparison."
            ),
            how_to_read=(
                f"The x-axis is {readout_channel} on a biexponential scale (linear near zero, log-like at high values); "
                f"the y-axis is % of each group's own peak count. Each line is one condition group, pooled across its "
                f"replicate wells. The dashed red line is the 1st-percentile-of-{control_label} gate; the dashed purple "
                f"line is the 50th-percentile-of-{control_label} gate."
            ),
        ))
        blocks.append(fig_hist)
        blocks.append(plot_overview(
            plot_start + 1, f"{title}: {readout_channel} vs SSC-A",
            why="Checks the same group separation in 2D, confirming the histogram shift isn't an artifact of pooling cells of very different granularity (SSC-A).",
            how="Same gated/compensated/pooled data as the histogram above, shown as a density contour (one contour per group) against SSC-A.",
            how_to_read="Each contour is one condition group's 2D density. The SSC-A axis is capped at the 99.5th percentile of plotted events (display only). Dashed lines mark the same two gates as the histogram.",
        ))
        blocks.append(fig_scatter)
        if len(df):
            blocks.append(plot_overview(
                plot_start + 2, f"{title}: % below 1st-percentile-of-mock gate",
                why=f"States the tail-effect knockdown number for {title} directly, per replicate well, with reader-friendly group/replicate labels.",
                how="Per-well % of gated cells falling below the 1st-percentile-of-mock gate.",
                how_to_read="Each bar is one replicate well, labeled by condition group and replicate number (never a well ID). Bar color matches the group.",
            ))
            blocks.append(bar_fig)
            blocks.append(plot_result(
                f"Mean knockdown for {primary_label} vs. {control_label} is {_kd1:+.1f} percentage points "
                f"(n={_n_primary:,} {primary_label} cells)."
            ))
            blocks.append(plot_overview(
                plot_start + 3, f"{title}: % below 50th-percentile-of-mock gate",
                why=f"States the population-level-shift knockdown number for {title} directly, per replicate well.",
                how="Per-well % of gated cells falling below the 50th-percentile-of-mock gate.",
                how_to_read="Same layout as the 1st-percentile bar chart above.",
            ))
            blocks.append(bar_fig_50)
            blocks.append(plot_result(
                f"Mean knockdown for {primary_label} vs. {control_label} is {_kd50:+.1f} percentage points "
                f"(n={_n_primary:,} {primary_label} cells)."
            ))
        if excl_note:
            blocks.append(mo.md(excl_note))
        if len(df):
            blocks.append(tidy_summary_table(df.sort_values(["label", "well"])))

        _cell_line = sample_sheet.loc[primary_wells[0], "Cell Line"] if primary_wells else ""
        summary = {
            "arm": title,
            "cell_line": _cell_line,
            "mock_group": control_label,
            "guide_group": primary_label,
            "knockdown_pp_1st": _kd1,
            "knockdown_pp_50th": _kd50,
            "flag": summary_flag,
        }
        return mo.vstack(blocks), summary


    return (build_arm,)


@app.cell(hide_code=True)
def _(build_arm, infection_gate_widgets):
    arm1_content_base, arm1_summary = build_arm(
        title="pDRT103/pDRT106 split-GFP, 8_3 cells",
        description=(
            "**Construct:** `pDRT103 <ORK|CD81>` is an all-in-one virus carrying both "
            "dCas9-KRAB and the sgRNA (ORK = non-targeting, CD81 = CD81-targeting); it "
            "requires co-infection with `pDRT106` to reconstitute split-GFP, so "
            "**GFP+ = confirmed dual infection**. Guide identity is read from which "
            "pDRT103 variant was used (B3/B4 = ORK, B5/B6 = CD81), not from a separate "
            "marker. Readout: CD81-BV421.\n\n"
            "**Three groups:** **ORK** (dCas9-KRAB + non-targeting guide -- transduction "
            "+ effector + guide control); **CD81** (dCas9-KRAB + CD81-targeting guide -- "
            "the real knockdown condition); and **pDRT109** (wells B1/B2, Sample IDs "
            "19/20) -- a transduction-only control vector that still requires `pDRT106` "
            "co-infection for the same split-GFP reconstitution (same GFP+ infection "
            "gate applies, no new marker needed), but carries **neither a guide nor "
            "dCas9-KRAB**. It controls for the physical act of viral transduction/"
            "integration itself, separate from ORK's guide+effector control -- expected "
            "to look mock-like (near-zero knockdown vs. ORK) if transduction alone "
            "doesn't perturb CD81 expression."
        ),
        groups=[
            ("ORK (pDRT103 ORK + pDRT106)", ["B3", "B4"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
            ("pDRT109 (transduction-only control + pDRT106)", ["B1", "B2"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
            ("CD81 (pDRT103 CD81 + pDRT106)", ["B5", "B6"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
        ],
        readout_channel="BV421-A",
        naive_well="A3",
        summary_flag="Includes pDRT109 transduction-only control (B1/B2, ~mock-like as expected) as a 3rd group in the panel below; not shown in this top-level ORK-vs-CD81 summary.",
        plot_start=19,
    )
    arm1_content_base
    return arm1_content_base, arm1_summary


@app.cell(hide_code=True)
def _(build_arm, infection_gate_widgets):
    arm2_content_base, arm2_summary = build_arm(
        title="AA173+AA239 two-virus, 17_3 cells",
        description=(
            "**Construct:** `AA173 BFP <ORK|CD81>` (guide-only, BFP-marked) delivered "
            "together with `AA239 dCas9-KRAB Thy1.1` (effector-only, Thy1.1-marked; "
            "detected here via the Thy1.1-FITC antibody stain since no FITC-tagged "
            "antibody is otherwise used in these wells). Functional knockdown requires "
            "**both** markers positive (BFP+ AND Thy1.1-FITC+). Readout: CD81-APC."
        ),
        groups=[
            ("ORK (AA173 ORK + AA239)", ["C5", "C6"],
             [("BV421-A", infection_gate_widgets.value["BV421-A__A2"]), ("FITC-A", infection_gate_widgets.value["FITC-A__A8_fmo"])]),
            ("CD81 (AA173 CD81 + AA239)", ["C7", "C8"],
             [("BV421-A", infection_gate_widgets.value["BV421-A__A2"]), ("FITC-A", infection_gate_widgets.value["FITC-A__A8_fmo"])]),
        ],
        readout_channel="APC-A",
        naive_well="A2",
        plot_start=23,
    )
    arm2_content_base
    return arm2_content_base, arm2_summary


@app.cell(hide_code=True)
def _(build_arm, infection_gate_widgets):
    arm3_content_base, arm3_summary = build_arm(
        title="AA173+AA239 two-virus, 8_3 cells",
        description=(
            "Same two-virus construct as the 17_3 panel above (BFP+ guide-only AA173 "
            "AND Thy1.1-FITC+ effector-only AA239, double-positive gate), but in 8_3 "
            "cells -- kept as its own panel rather than pooled with the 17_3 arm since "
            "cell line is a separate biological background. Readout: CD81-APC."
        ),
        groups=[
            ("ORK (AA173 ORK + AA239)", ["C9", "C10"],
             [("BV421-A", infection_gate_widgets.value["BV421-A__A3"]), ("FITC-A", infection_gate_widgets.value["FITC-A__D5_fmo"])]),
            ("CD81 (AA173 CD81 + AA239)", ["C11", "C12"],
             [("BV421-A", infection_gate_widgets.value["BV421-A__A3"]), ("FITC-A", infection_gate_widgets.value["FITC-A__D5_fmo"])]),
        ],
        readout_channel="APC-A",
        naive_well="A3",
        plot_start=27,
    )
    arm3_content_base
    return arm3_content_base, arm3_summary


@app.cell(hide_code=True)
def _(build_arm, infection_gate_widgets):
    arm4_content_base, arm4_summary = build_arm(
        title="AA228 mCherry all-in-one, 8_3 cells",
        description=(
            "**Construct:** `AA228 mCherry <ORK|CD81>` is an all-in-one virus (guide + "
            "dCas9-KRAB together), marked with mCherry -- single-marker infection gate "
            "(mCherry+, `Y610-mCHERRY-A` channel). Readout: CD81-BV421."
        ),
        groups=[
            ("ORK (AA228 mCherry ORK)", ["G4", "G5"], [("Y610-mCHERRY-A", infection_gate_widgets.value["Y610-mCHERRY-A__A3"])]),
            ("CD81 (AA228 mCherry CD81)", ["G6", "G7"], [("Y610-mCHERRY-A", infection_gate_widgets.value["Y610-mCHERRY-A__A3"])]),
        ],
        readout_channel="BV421-A",
        naive_well="A3",
        summary_flag="Replicate split: G7 shows strong tail effect, G6 does not (see bar charts).",
        plot_start=31,
    )
    arm4_content_base
    return arm4_content_base, arm4_summary


@app.cell(hide_code=True)
def _(build_arm, infection_gate_widgets):
    arm5_content, arm5_summary = build_arm(
        title="WTC11-KRAB integrated line + AA173 guide, Dox 2x2",
        description=(
            "**Construct:** dCas9-KRAB is stably integrated and Dox-inducible in this "
            "WTC11 sub-line; the guide is delivered lentivirally via `AA173 BFP "
            "<ORK|CD81>` (BFP+ = guide-delivery gate). A clean 2x2: guide identity "
            "(ORK: H3/H4, CD81: H5-H8) x Dox induction (No: H3,H5,H6 / Yes: H4,H7,H8). "
            "**Primary comparison** is Dox-on (H7,H8) vs Dox-off (H5,H6) *within* the "
            "CD81 guide; ORK wells (H3 vs H4) are a secondary guide-identity reference "
            "at each Dox state. Readout: CD81-APC. Gate is fixed at the median of the "
            "Dox-off ORK group (H3) so all four groups are compared on the same scale."
        ),
        groups=[
            ("ORK, Dox-off", ["H3"], [("BV421-A", infection_gate_widgets.value["BV421-A__A1"])]),
            ("ORK, Dox-on", ["H4"], [("BV421-A", infection_gate_widgets.value["BV421-A__A1"])]),
            ("CD81, Dox-off", ["H5", "H6"], [("BV421-A", infection_gate_widgets.value["BV421-A__A1"])]),
            ("CD81, Dox-on", ["H7", "H8"], [("BV421-A", infection_gate_widgets.value["BV421-A__A1"])]),
        ],
        readout_channel="APC-A",
        naive_well="A1",
        control_group_index=0,
        plot_start=35,
    )
    arm5_content
    return arm5_content, arm5_summary


@app.cell(hide_code=True)
def _(build_arm, infection_gate_widgets):
    arm6_content_base, arm6_summary = build_arm(
        title="WTC11 fully-integrated CD81-BFP guide+effector, C1-C4 only",
        description=(
            "All four wells (C1-C4) share the same WTC11 background with the "
            "CD81-targeting guide stably integrated (`Integrated Guide = CD81-BFP`, "
            "BFP+ = integrated-guide gate) and no Dox induction. The comparison is "
            "**effector delivery**, not guide identity: C3/C4 received no additional "
            "virus (no effector on top of the integrated guide), while C1/C2 also "
            "received `AA239 dCas9-KRAB Thy1.1` (Thy1.1-FITC+ = confirmed effector "
            "delivery). Checked directly: C3/C4 show only ~2-4% Thy1.1+ (background "
            "level, consistent with no AA239 virus delivered), so C3/C4 are gated on "
            "BFP+ only -- a Thy1.1+ gate would be meaningless (and would discard "
            "nearly all cells) for wells that never received the Thy1.1-marked virus. "
            "C1/C2 are gated on the full double-positive BFP+ AND Thy1.1+ population, "
            "same idiom as the AA173+AA239 arms above. The 1st-percentile-of-mock gate "
            "is set from the no-effector group (C3/C4). Readout: CD81-APC."
        ),
        groups=[
            ("No effector (BFP+ only, C3/C4)", ["C3", "C4"], [("BV421-A", infection_gate_widgets.value["BV421-A__A1"])]),
            ("Effector delivered (BFP+ AND Thy1.1+, C1/C2)", ["C1", "C2"],
             [("BV421-A", infection_gate_widgets.value["BV421-A__A1"]), ("FITC-A", infection_gate_widgets.value["FITC-A__A4_fmo"])]),
        ],
        readout_channel="APC-A",
        naive_well="A1",
        control_group_index=0,
        summary_flag="Redefined: C3/C4 (no effector) vs C1/C2 (+AA239 effector), not the original Dox=Yes/No design.",
        plot_start=39,
    )
    arm6_content_base
    return arm6_content_base, arm6_summary


@app.cell
def _(mo):
    ssc_cap_slider = mo.ui.slider(
        start=0, stop=5_000_000, step=10_000, value=3_000_000,
        label="SSC-A upper cap (cell gate; events above this SSC-A value are excluded)",
        full_width=True, show_value=True,
    )
    ssc_cap_slider
    return (ssc_cap_slider,)


@app.cell
def _(WELLS, debris_gate, go, mo, pd, plot_overview, plot_result):
    _TIME_QC_ROWS = []
    for _w in WELLS:
        _d = debris_gate(_w)
        if len(_d) < 50:
            continue
        _bins = pd.qcut(_d["Time"], q=10, duplicates="drop")
        _medians = _d.groupby(_bins, observed=True)["FSC-A"].median()
        _swing = (float(_medians.max()) - float(_medians.min())) / max(float(_medians.median()), 1.0)
        _TIME_QC_ROWS.append({"Well": _w, "n (gated)": len(_d), "FSC-A median swing across time bins": round(100 * _swing, 1)})

    _time_qc_df = pd.DataFrame(_TIME_QC_ROWS).sort_values("FSC-A median swing across time bins", ascending=False)
    _flagged = _time_qc_df[_time_qc_df["FSC-A median swing across time bins"] > 25]

    _time_fig = go.Figure()
    _time_fig.add_trace(go.Bar(
        x=_time_qc_df["Well"], y=_time_qc_df["FSC-A median swing across time bins"],
        marker_color=["#E45756" if v > 25 else "#4C78A8" for v in _time_qc_df["FSC-A median swing across time bins"]],
    ))
    _time_fig.add_shape(type="line", x0=0, x1=1, xref="paper", y0=25, y1=25, line=dict(color="red", width=2, dash="dash"))
    _time_fig.update_layout(
        title="Plot 18. Acquisition-stability check: FSC-A median swing across 10 time bins, per well",
        yaxis_title="% swing (max-min median FSC-A / overall median)", height=420, margin=dict(t=60),
    )

    time_qc_section = mo.vstack([
        plot_overview(
            18, "Acquisition-stability check",
            why="Checks for acquisition instability (a clog or bubble partway through a well's run) that would silently bias that well's gated population, independent of any biological effect.",
            how="Each well's Time channel is split into 10 equal-count bins; the median FSC-A is computed per bin, and the swing is (max bin median - min bin median) / overall median FSC-A.",
            how_to_read="One bar per well; the dashed red line marks a 25% swing. Bars above the line are flagged below as a caution, not excluded by a hard filter unless a well shows a sharp jump rather than a mild swing.",
        ),
        _time_fig,
        plot_result(
            (
                f"{len(_flagged)} of {len(_time_qc_df)} wells exceed a 25% FSC-A median swing across time bins: "
                + ", ".join(
                    f"{row['Well']} ({row['FSC-A median swing across time bins']:.0f}%)"
                    for _, row in _flagged.iterrows()
                )
            ) if len(_flagged) else
            f"0 of {len(_time_qc_df)} wells exceed a 25% FSC-A median swing across time bins -- no acquisition-stability concerns flagged."
        ),
    ])
    time_qc_section
    return


if __name__ == "__main__":
    app.run()
