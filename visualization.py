from pathlib import Path

import matplotlib
import numpy as np


matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
from plotly.subplots import make_subplots


DATA_PATH = Path("dataRaw.csv")
OUTPUT_DIR = Path("outputs")
PLOTLY_OUTPUT_DIR = Path("outputs_plotly")
VELOCITIES = [3.5, 3.0, 2.5]
CONCENTRATIONS = [0.1, 0.2]
EARLY_WINDOW_MIN = 10.0
LATE_WINDOW_MIN = 10.0
TEMPERATURE_REF_C = 25.0
CURRENT_COMPENSATION_PER_C = 0.01
CONDUCTIVITY_COMPENSATION_PER_C = 0.02
TEMPERATURE_REFERENCE_COL = "TIT1"
CURRENT_25_COL = "I25"
RESISTANCE_25_COL = "R25"
POWER_25_COL = "P25"
TIME_TICKS_MIN = np.arange(0, 91, 15)
OBSOLETE_OUTPUTS = [
    "delta_rcalc_profiles.png",
    "limiting_proxy_profiles.png",
    "limiting_proxy.csv",
    "resistance_ohm_article_style.png",
    "resistance_vs_inverse_current.png",
    "resistance_calc_article_style.png",
    "resistance_calc_vs_inverse_current.png",
    "run_30_45_current_resistance_overview.png",
    "run_30_45_diluate_concentrate_overview.png",
]
RUN_SELECTIONS = [(1, 14), (29, 45)]
FULL_RUN_SELECTIONS = [(1, 45)]
HYDRAULICS_AVG_WINDOW_MIN = 5
PLOTLY_RUN_IDS = [15, 17, 19, 21]
PLOTLY_CONDUCTIVITY_COLS = ["CO1_25", "CO2_25", "CO3_25", "CO4_25", "CO5_25"]


def load_data(path: Path) -> pd.DataFrame:
    """Load the raw CSV and create temperature-compensated working columns."""
    df = pd.read_csv(path)
    df["date/time"] = pd.to_datetime(df["date/time"], errors="coerce")
    current_factor = 1.0 + CURRENT_COMPENSATION_PER_C * (df[TEMPERATURE_REFERENCE_COL] - TEMPERATURE_REF_C)
    conductivity_factor = 1.0 + CONDUCTIVITY_COMPENSATION_PER_C * (df[TEMPERATURE_REFERENCE_COL] - TEMPERATURE_REF_C)
    df[CURRENT_25_COL] = np.where(current_factor.abs() > 1e-9, df["I"] / current_factor, np.nan)
    df[RESISTANCE_25_COL] = np.where(df[CURRENT_25_COL].abs() > 1e-9, df["Ustack"] / df[CURRENT_25_COL], np.nan)
    df[POWER_25_COL] = df["Ustack"] * df[CURRENT_25_COL]
    for column in ["CO1", "CO2", "CO3", "CO4", "CO5"]:
        df[f"{column}_25"] = np.where(conductivity_factor.abs() > 1e-9, df[column] / conductivity_factor, np.nan)
    return df


def cleanup_outputs(output_dir: Path) -> None:
    """Remove legacy generated files that are no longer part of the workflow."""
    for filename in OBSOLETE_OUTPUTS:
        path = output_dir / filename
        if path.exists():
            path.unlink()


def add_run_id(df: pd.DataFrame) -> pd.DataFrame:
    """Assign a run identifier whenever time, flow condition, concentration, or current resets."""
    run_reset = (
        df["Time per step (min)"].diff().fillna(0).lt(0)
        | df["v"].diff().fillna(0).ne(0)
        | df["Ce"].diff().fillna(0).ne(0)
        | df["I"].diff().fillna(0).ne(0)
    )
    enriched = df.copy()
    enriched["run_id"] = run_reset.cumsum() + 1
    return enriched


def prepare_active_runs(df: pd.DataFrame) -> pd.DataFrame:
    """Filter to positive-current runs at the target velocities and add display labels."""
    active = df.loc[(df["I"] > 0) & (df["v"].isin(VELOCITIES))].copy()
    active["v_label"] = active["v"].map(lambda value: f"{value:.1f} cm/s")
    active["Ce_label"] = active["Ce"].map(lambda value: f"{value:.1f} M")
    active["I_label"] = active["I"].map(lambda value: f"{value:.2f} A")
    return active


def _selection_label(selections: list[tuple[int, int]]) -> str:
    """Convert run-range selections to a filename-friendly label."""
    return "_".join(f"{start}_{end}" for start, end in selections)


def _selection_title(selections: list[tuple[int, int]]) -> str:
    """Convert run-range selections to a human-readable title fragment."""
    return ", ".join(f"{start}-{end}" for start, end in selections)


def prepare_run_window(df: pd.DataFrame, selections: list[tuple[int, int]]) -> pd.DataFrame:
    """Build a stitched cumulative-time block for the selected run ranges."""
    mask = pd.Series(False, index=df.index)
    for start_run, end_run in selections:
        mask |= df["run_id"].between(start_run, end_run)

    window = df.loc[mask].copy()
    if window.empty:
        return window

    window = window.sort_values(["run_id", "Time per step (min)", "date/time"]).copy()
    cumulative_offsets_min: dict[int, float] = {}
    offset_min = 0.0
    run_order = sorted(window["run_id"].unique().tolist())
    for run_id in run_order:
        run_slice = window.loc[window["run_id"] == run_id, "Time per step (min)"]
        cumulative_offsets_min[run_id] = offset_min
        offset_min += run_slice.max()

    window["cumulative_time_h"] = (
        window["run_id"].map(cumulative_offsets_min) + window["Time per step (min)"]
    ) / 60.0
    window["v_label"] = window["v"].map(lambda value: f"{value:.1f} cm/s")
    window["delta_rlinear"] = window["Rlinear"] - window.groupby("run_id")["Rlinear"].transform("first")
    return window


def average_run_window_signals(
    run_window: pd.DataFrame,
    signal_columns: list[str],
    window_min: int,
) -> pd.DataFrame:
    """Average selected signals into fixed-width time bins within each run."""
    averaged = run_window.copy()
    averaged["time_bin_min"] = (
        np.floor(averaged["cumulative_time_h"] * 60.0 / window_min) * window_min
    )
    grouped = averaged.groupby(["run_id", "v", "v_label", "time_bin_min"], as_index=False)[signal_columns].mean()
    grouped["cumulative_time_h"] = (grouped["time_bin_min"] + window_min / 2.0) / 60.0
    return grouped


def add_relative_metrics(active: pd.DataFrame) -> pd.DataFrame:
    """Add run-relative metrics used by the comparison plots."""
    enriched = active.copy()
    grouped = enriched.groupby("run_id")
    enriched["delta_rlinear"] = enriched["Rlinear"] - grouped["Rlinear"].transform("first")
    enriched["norm_delta_rlinear"] = enriched["delta_rlinear"] / grouped["Rlinear"].transform("first")
    enriched["inv_I25"] = 1.0 / enriched[CURRENT_25_COL]
    return enriched


def prepare_plotly_runs(df: pd.DataFrame, run_ids: list[int]) -> pd.DataFrame:
    """Prepare the selected runs for standalone Plotly HTML pages."""
    selected = df.loc[df["run_id"].isin(run_ids)].copy()
    if selected.empty:
        return selected

    grouped = selected.groupby("run_id")
    selected["delta_rlinear"] = selected["Rlinear"] - grouped["Rlinear"].transform("first")
    selected["run_label"] = selected["run_id"].map(lambda value: f"Run {value}")
    return selected.sort_values(["run_id", "Time per step (min)", "date/time"]).copy()


def build_run_summary(active: pd.DataFrame) -> pd.DataFrame:
    """Summarize each active run with start/end resistance and operating metadata."""
    grouped = active.sort_values("Time per step (min)").groupby("run_id", as_index=False)

    summary = grouped.agg(
        v=("v", "first"),
        Ce=("Ce", "first"),
        I=("I", "first"),
        I25=(CURRENT_25_COL, "first"),
        duration_min=("Time per step (min)", "max"),
        final_rlinear=("Rlinear", "last"),
        initial_rlinear=("Rlinear", "first"),
        final_r25=(RESISTANCE_25_COL, "last"),
        initial_r25=(RESISTANCE_25_COL, "first"),
        final_ustack=("Ustack", "last"),
    )
    summary["delta_rlinear"] = summary["final_rlinear"] - summary["initial_rlinear"]
    summary["v_label"] = summary["v"].map(lambda value: f"{value:.1f} cm/s")
    summary["Ce_label"] = summary["Ce"].map(lambda value: f"{value:.1f} M")
    return summary


def build_all_run_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Summarize every run and classify it by polarity and duration."""
    summary = (
        df.sort_values("Time per step (min)")
        .groupby("run_id", as_index=False)
        .agg(
            v=("v", "first"),
            Ce=("Ce", "first"),
            I=("I", "first"),
            duration_min=("Time per step (min)", "max"),
        )
    )
    summary["status"] = np.where(
        summary["I"] < 0,
        "negative",
        np.where(summary["duration_min"] < 50, "active_short", "active"),
    )
    return summary


def build_window_summary(active: pd.DataFrame) -> pd.DataFrame:
    """Compute early and late window medians for voltage and resistance metrics."""
    rows = []
    for run_id, run_df in active.groupby("run_id"):
        run_df = run_df.sort_values("Time per step (min)")
        max_time = run_df["Time per step (min)"].max()
        early_mask = run_df["Time per step (min)"] <= EARLY_WINDOW_MIN
        late_mask = run_df["Time per step (min)"] >= max(0.0, max_time - LATE_WINDOW_MIN)
        rows.append(
            {
                "run_id": run_id,
                "v": run_df["v"].iat[0],
                "Ce": run_df["Ce"].iat[0],
                "I": run_df["I"].iat[0],
                "duration_min": max_time,
                "u_early": run_df.loc[early_mask, "Ustack"].median(),
                "u_late": run_df.loc[late_mask, "Ustack"].median(),
                "r_early": run_df.loc[early_mask, "Rlinear"].median(),
                "r_late": run_df.loc[late_mask, "Rlinear"].median(),
                "r25_early": run_df.loc[early_mask, RESISTANCE_25_COL].median(),
                "r25_late": run_df.loc[late_mask, RESISTANCE_25_COL].median(),
            }
        )

    summary = pd.DataFrame(rows)
    summary["v_label"] = summary["v"].map(lambda value: f"{value:.1f} cm/s")
    summary["Ce_label"] = summary["Ce"].map(lambda value: f"{value:.1f} M")
    return summary


def build_inverse_current_summary(active: pd.DataFrame, summary: pd.DataFrame, window_summary: pd.DataFrame) -> pd.DataFrame:
    """Prepare phase-based run summaries for reciprocal-current plots."""
    resistance_summary = (
        active.sort_values("Time per step (min)")
        .groupby("run_id", as_index=False)
        .agg(
            initial_r25=(RESISTANCE_25_COL, "first"),
            final_r25=(RESISTANCE_25_COL, "last"),
            initial_rlinear=("Rlinear", "first"),
            final_rlinear=("Rlinear", "last"),
            initial_i25=(CURRENT_25_COL, "first"),
        )
    )

    merged = summary.merge(
        window_summary[["run_id", "u_early", "u_late", "r_early", "r_late", "r25_early", "r25_late"]],
        on="run_id",
        how="left",
    ).merge(resistance_summary, on="run_id", how="left", suffixes=("", "_dup"))
    merged["inv_I25"] = 1.0 / merged["I25"]

    rows = []
    phase_specs = [
        ("start", "initial_r25", "initial_rlinear"),
        ("early", "r25_early", "r_early"),
        ("late", "r25_late", "r_late"),
        ("end", "final_r25", "final_rlinear"),
    ]
    for _, row in merged.iterrows():
        for phase, resistance_col, rlinear_col in phase_specs:
            rows.append(
                {
                    "run_id": row["run_id"],
                    "v": row["v"],
                    "Ce": row["Ce"],
                    "I": row["I"],
                    "I25": row["I25"],
                    "inv_I25": row["inv_I25"],
                    "phase": phase,
                    "phase_order": ["start", "early", "late", "end"].index(phase),
                    "r25": row[resistance_col],
                    "rlinear": row[rlinear_col],
                    "v_label": row["v_label"],
                    "Ce_label": row["Ce_label"],
                }
            )

    return pd.DataFrame(rows)


def apply_time_axis(grid: sns.FacetGrid) -> None:
    """Apply a consistent 15-minute x-axis to time-domain facet plots."""
    for ax in grid.axes.flat:
        ax.set_xticks(TIME_TICKS_MIN)
        ax.set_xlim(TIME_TICKS_MIN[0], TIME_TICKS_MIN[-1])


def annotate_missing_facets(
    grid: sns.FacetGrid,
    data: pd.DataFrame,
    row_key: str = "v_label",
    col_key: str = "Ce_label",
) -> None:
    """Mark empty facets for untested combinations such as missing velocity/concentration pairs."""
    axes = np.array(grid.axes, dtype=object).reshape(len(grid.row_names), len(grid.col_names))
    for row_index, row_value in enumerate(grid.row_names):
        for col_index, col_value in enumerate(grid.col_names):
            ax = axes[row_index, col_index]
            facet_data = data.loc[(data[row_key] == row_value) & (data[col_key] == col_value)]
            if facet_data.empty:
                ax.text(
                    0.5,
                    0.5,
                    "Not tested",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=12,
                    color="dimgray",
                )


def save_profiles(active: pd.DataFrame) -> None:
    """Save the main faceted time-series profiles for voltage and resistance metrics."""
    profile_specs = [
        ("Ustack", "Stack voltage [V]", "voltage_article_style.png"),
        (RESISTANCE_25_COL, "R25 [Ohm]", "resistance_25_article_style.png"),
        ("Rlinear", "Linearized resistance", "rlinear_profiles.png"),
        ("Ustack", "Stack voltage [V]", "ustack_profiles.png"),
        ("delta_rlinear", "Delta Rlinear", "delta_rlinear_profiles.png"),
        ("norm_delta_rlinear", "Relative delta Rlinear", "norm_delta_rlinear_profiles.png"),
    ]

    for y_col, y_label, filename in profile_specs:
        grid = sns.relplot(
            data=active.sort_values("Time per step (min)"),
            x="Time per step (min)",
            y=y_col,
            hue="I_label",
            kind="line",
            units="run_id",
            estimator=None,
            row="v_label",
            col="Ce_label",
            row_order=[f"{value:.1f} cm/s" for value in VELOCITIES],
            col_order=[f"{value:.1f} M" for value in CONCENTRATIONS],
            height=3.2,
            aspect=1.2,
            linewidth=2.0,
            palette="viridis",
            facet_kws={"sharey": False, "margin_titles": True},
        )
        grid.set_axis_labels("time [min]", y_label)
        apply_time_axis(grid)
        grid.figure.subplots_adjust(top=0.9)
        if filename == "resistance_25_article_style.png":
            grid.figure.suptitle("R25 [Ohm] = f(time) [min]")
        else:
            grid.figure.suptitle(f"{y_label} in function of time")
        annotate_missing_facets(grid, active)
        grid.savefig(OUTPUT_DIR / filename, dpi=200, bbox_inches="tight")
        plt.close(grid.figure)


def save_run_summary(summary: pd.DataFrame) -> None:
    """Save the run-level scatter summary used for quick condition comparison."""
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    sns.scatterplot(
        data=summary,
        x="I",
        y="delta_rlinear",
        hue="v_label",
        style="Ce_label",
        size="final_ustack",
        sizes=(70, 260),
        palette="crest",
        ax=ax,
    )
    ax.set_xlabel("Applied current [A]")
    ax.set_ylabel("Delta Rlinear over active cycle")
    ax.set_title("Run-level fouling summary")
    ax.grid(alpha=0.25)
    if ax.legend_ is not None:
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0.0, frameon=True)
    fig.savefig(OUTPUT_DIR / "run_summary_scatter.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_run_duration_status_summary(all_run_summary: pd.DataFrame) -> None:
    """Save a status overview of run durations across active and negative-polarity runs."""
    fig, ax = plt.subplots(figsize=(10, 5.5))

    negative = all_run_summary.loc[all_run_summary["status"] == "negative"]
    active = all_run_summary.loc[all_run_summary["status"] == "active"]
    active_short = all_run_summary.loc[all_run_summary["status"] == "active_short"]

    if not negative.empty:
        ax.scatter(
            negative["run_id"],
            negative["duration_min"],
            color="royalblue",
            marker="_",
            s=220,
            linewidths=2.2,
            label="Negative polarity",
        )
    if not active.empty:
        ax.scatter(
            active["run_id"],
            active["duration_min"],
            color="red",
            marker="+",
            s=180,
            linewidths=2.2,
            label="Active >= 50 min",
        )
    if not active_short.empty:
        ax.scatter(
            active_short["run_id"],
            active_short["duration_min"],
            color="orange",
            marker="+",
            s=180,
            linewidths=2.2,
            label="Active < 50 min",
        )

    ax.set_xlabel("Run ID")
    ax.set_ylabel("Duration [min]")
    ax.set_title("Run duration by run ID")
    ax.grid(alpha=0.25)
    run_ticks = [1] + list(range(5, int(all_run_summary["run_id"].max()) + 1, 5))
    ax.set_xticks(run_ticks)
    if ax.legend_ is not None:
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0.0, frameon=True)
    fig.savefig(OUTPUT_DIR / "run_duration_status_scatter.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_inverse_current_plots(inverse_current: pd.DataFrame) -> None:
    """Save the reciprocal-current plot for Rlinear under the 0.1 M condition."""
    phase_order = ["start", "early", "late", "end"]
    rlinear_data = inverse_current.loc[inverse_current["Ce"] == 0.1].copy()
    if not rlinear_data.empty:
        grid = sns.relplot(
            data=rlinear_data.sort_values(["phase_order", "inv_I25"]),
            x="inv_I25",
            y="rlinear",
            hue="phase",
            style="phase",
            hue_order=phase_order,
            style_order=phase_order,
            kind="line",
            markers=True,
            dashes=False,
            row="v_label",
            col="Ce_label",
            row_order=[f"{value:.1f} cm/s" for value in VELOCITIES],
            col_order=["0.1 M"],
            height=3.2,
            aspect=1.2,
            facet_kws={"sharey": False, "margin_titles": True},
        )
        grid.set_axis_labels("1 / I25 [A$^{-1}$]", "Rlinear")
        grid.figure.subplots_adjust(top=0.9)
        grid.figure.suptitle("Rlinear vs reciprocal compensated current")
        annotate_missing_facets(grid, rlinear_data)
        grid.savefig(OUTPUT_DIR / "rlinear_vs_inverse_current.png", dpi=200, bbox_inches="tight")
        plt.close(grid.figure)


def save_resistance_25_inverse_current(active: pd.DataFrame) -> None:
    """Save the all-points reciprocal-current plot for temperature-normalized resistance."""
    grid = sns.relplot(
        data=active.sort_values(["v", "Ce", "inv_I25", "Time per step (min)"]),
        x="inv_I25",
        y=RESISTANCE_25_COL,
        hue="Time per step (min)",
        kind="scatter",
        row="v_label",
        col="Ce_label",
        row_order=[f"{value:.1f} cm/s" for value in VELOCITIES],
        col_order=[f"{value:.1f} M" for value in CONCENTRATIONS],
        height=3.2,
        aspect=1.2,
        palette="viridis",
        alpha=0.35,
        s=18,
        legend=False,
        facet_kws={"sharey": False, "margin_titles": True},
    )
    grid.set_axis_labels("1 / I25 [A$^{-1}$]", "R25 [Ohm]")
    grid.figure.subplots_adjust(top=0.9)
    grid.figure.suptitle("R25 [Ohm] vs reciprocal compensated current")
    annotate_missing_facets(grid, active)
    grid.savefig(OUTPUT_DIR / "resistance_25_vs_inverse_current.png", dpi=200, bbox_inches="tight")
    plt.close(grid.figure)


def save_run_window_overview(run_window: pd.DataFrame, selections: list[tuple[int, int]]) -> None:
    """Save the 2x2 stitched overview for a selected run block."""
    if run_window.empty:
        return

    plot_data = run_window.copy()
    selection_title = _selection_title(selections)
    selection_label = _selection_label(selections)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    palette = sns.color_palette("viridis", n_colors=run_window["v_label"].nunique())

    sns.lineplot(
        data=plot_data.sort_values(["run_id", "cumulative_time_h"]),
        x="cumulative_time_h",
        y=CURRENT_25_COL,
        hue="v_label",
        units="run_id",
        estimator=None,
        palette=palette,
        linewidth=1.8,
        ax=axes[0, 0],
    )
    axes[0, 0].set_xlabel("Cumulative time [h]")
    axes[0, 0].set_ylabel("I25 [A]")
    axes[0, 0].set_title(f"Runs {selection_title}: compensated current")
    axes[0, 0].grid(alpha=0.25)

    sns.lineplot(
        data=plot_data.sort_values(["run_id", "cumulative_time_h"]),
        x="cumulative_time_h",
        y=RESISTANCE_25_COL,
        hue="v_label",
        units="run_id",
        estimator=None,
        palette=palette,
        linewidth=1.8,
        legend=False,
        ax=axes[0, 1],
    )
    axes[0, 1].set_xlabel("Cumulative time [h]")
    axes[0, 1].set_ylabel("R25 [Ohm]")
    axes[0, 1].set_title(f"Runs {selection_title}: compensated resistance")
    axes[0, 1].grid(alpha=0.25)

    sns.lineplot(
        data=plot_data.sort_values(["run_id", "cumulative_time_h"]),
        x="cumulative_time_h",
        y="TIT1",
        hue="v_label",
        units="run_id",
        estimator=None,
        palette=palette,
        linewidth=1.8,
        legend=False,
        ax=axes[1, 0],
    )
    axes[1, 0].set_xlabel("Cumulative time [h]")
    axes[1, 0].set_ylabel("TIT1 [C]")
    axes[1, 0].set_title(f"Runs {selection_title}: temperature TIT1")
    axes[1, 0].grid(alpha=0.25)

    sns.lineplot(
        data=plot_data.sort_values(["run_id", "cumulative_time_h"]),
        x="cumulative_time_h",
        y=POWER_25_COL,
        hue="v_label",
        units="run_id",
        estimator=None,
        palette=palette,
        linewidth=1.8,
        legend=False,
        ax=axes[1, 1],
    )
    axes[1, 1].set_xlabel("Cumulative time [h]")
    axes[1, 1].set_ylabel("P25 [W]")
    axes[1, 1].set_title(f"Runs {selection_title}: compensated power")
    axes[1, 1].grid(alpha=0.25)

    if axes[0, 0].legend_ is not None:
        sns.move_legend(axes[0, 0], "upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0.0, frameon=True)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f"run_{selection_label}_current_resistance_overview.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_run_window_hydraulics_overview(run_window: pd.DataFrame, selections: list[tuple[int, int]]) -> None:
    """Save the 3x2 stitched hydraulics and resistance overview for a selected run block."""
    if run_window.empty:
        return

    plot_data = run_window.sort_values(["run_id", "cumulative_time_h"]).copy()
    selection_title = _selection_title(selections)
    selection_label = _selection_label(selections)
    hydraulics_data = average_run_window_signals(
        plot_data,
        signal_columns=["FIT2", "FIT1", "PIT4", "PIT3"],
        window_min=HYDRAULICS_AVG_WINDOW_MIN,
    )
    fig, axes = plt.subplots(3, 2, figsize=(13, 12), sharex=True)
    palette = sns.color_palette("viridis", n_colors=plot_data["v_label"].nunique())

    panel_specs = [
        ("FIT2", "FIT2 [l/h]", f"Runs {selection_title}: diluate flow ({HYDRAULICS_AVG_WINDOW_MIN}-min mean)", hydraulics_data),
        ("FIT1", "FIT1 [l/h]", f"Runs {selection_title}: concentrate flow ({HYDRAULICS_AVG_WINDOW_MIN}-min mean)", hydraulics_data),
        ("PIT4", "PIT4 [bar]", f"Runs {selection_title}: diluate pressure ({HYDRAULICS_AVG_WINDOW_MIN}-min mean)", hydraulics_data),
        ("PIT3", "PIT3 [bar]", f"Runs {selection_title}: concentrate pressure ({HYDRAULICS_AVG_WINDOW_MIN}-min mean)", hydraulics_data),
        ("delta_rlinear", "Delta Rlinear", f"Runs {selection_title}: resistance increase metric", plot_data),
        ("TIT1", "TIT1 [C]", f"Runs {selection_title}: temperature", plot_data),
    ]

    for ax, (y_col, y_label, title, data_source) in zip(axes.flat, panel_specs):
        use_legend = y_col == "FIT2"
        sns.lineplot(
            data=data_source,
            x="cumulative_time_h",
            y=y_col,
            hue="v_label",
            units="run_id",
            estimator=None,
            palette=palette,
            linewidth=1.8,
            legend=use_legend,
            ax=ax,
        )
        ax.set_xlabel("Cumulative time [h]")
        ax.set_ylabel(y_label)
        ax.set_title(title)
        ax.grid(alpha=0.25)
        if use_legend and ax.legend_ is not None:
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0.0, frameon=True)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f"run_{selection_label}_diluate_concentrate_overview.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def export_tables(
    summary: pd.DataFrame,
    window_summary: pd.DataFrame,
    inverse_current: pd.DataFrame,
) -> None:
    """Write the summary tables consumed by the static analysis workflow."""
    summary.sort_values(["Ce", "v", "I"]).to_csv(OUTPUT_DIR / "run_summary.csv", index=False)
    window_summary.sort_values(["Ce", "v", "I"]).to_csv(OUTPUT_DIR / "window_summary.csv", index=False)
    inverse_current.sort_values(["Ce", "v", "phase_order", "I"]).to_csv(
        OUTPUT_DIR / "inverse_current_summary.csv",
        index=False,
    )


def _format_timestamp(value: pd.Timestamp) -> str:
    """Format timestamps for the Plotly index page."""
    if pd.isna(value):
        return "n/a"
    return value.strftime("%Y-%m-%d %H:%M:%S")


def _build_plotly_index(output_dir: Path, pages: list[dict[str, str]]) -> None:
    """Write a small HTML index for the generated Plotly run pages."""
    if not pages:
        return

    rows = "\n".join(
        (
            "<tr>"
            f"<td>{page['run_id']}</td>"
            f"<td>{page['v']}</td>"
            f"<td>{page['Ce']}</td>"
            f"<td>{page['I']}</td>"
            f"<td>{page['duration_min']}</td>"
            f"<td>{page['start_time']}</td>"
            f"<td>{page['end_time']}</td>"
            f"<td><a href=\"{page['filename']}\">{page['filename']}</a></td>"
            "</tr>"
        )
        for page in pages
    )
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Selected Run Plotly Exports</title>
  <style>
    body {{
      font-family: "Segoe UI", Arial, sans-serif;
      margin: 2rem;
      background: #f5f7fa;
      color: #1f2933;
    }}
    h1 {{
      margin-bottom: 0.5rem;
    }}
    p {{
      max-width: 70rem;
    }}
    table {{
      border-collapse: collapse;
      width: 100%;
      background: #ffffff;
      box-shadow: 0 8px 24px rgba(15, 23, 42, 0.08);
    }}
    th, td {{
      padding: 0.8rem;
      border-bottom: 1px solid #d9e2ec;
      text-align: left;
    }}
    th {{
      background: #e9eef5;
    }}
    a {{
      color: #0b69a3;
      text-decoration: none;
    }}
  </style>
</head>
<body>
  <h1>Plotly Exports For Selected Runs</h1>
  <p>Interactive HTML pages generated from <code>visualization.py</code> for runs 15, 17, 19, and 21.</p>
  <table>
    <thead>
      <tr>
        <th>Run ID</th>
        <th>Velocity [cm/s]</th>
        <th>Ce [M]</th>
        <th>I [A]</th>
        <th>Duration [min]</th>
        <th>Start time</th>
        <th>End time</th>
        <th>Page</th>
      </tr>
    </thead>
    <tbody>
      {rows}
    </tbody>
  </table>
</body>
</html>
"""
    (output_dir / "index.html").write_text(html, encoding="utf-8")


def save_plotly_run_pages(selected_runs: pd.DataFrame, output_dir: Path) -> None:
    """Save standalone Plotly HTML pages for the requested run IDs."""
    if selected_runs.empty:
        return

    output_dir.mkdir(exist_ok=True)
    page_rows: list[dict[str, str]] = []
    conductivity_colors = ["#0f766e", "#1d4ed8", "#9333ea", "#c2410c", "#ca8a04"]

    for run_id in PLOTLY_RUN_IDS:
        run_df = selected_runs.loc[selected_runs["run_id"] == run_id].copy()
        if run_df.empty:
            continue

        run_df = run_df.sort_values(["Time per step (min)", "date/time"]).copy()
        run_meta = run_df.iloc[0]
        filename = f"run_{run_id}_interactive.html"
        duration_min = run_df["Time per step (min)"].max()
        title = (
            f"Run {run_id} interactive overview"
            f"<br><sup>v={run_meta['v']:.1f} cm/s | Ce={run_meta['Ce']:.1f} M | "
            f"I={run_meta['I']:.2f} A | duration={duration_min:.2f} min</sup>"
        )

        fig = make_subplots(
            rows=3,
            cols=2,
            specs=[
                [{"secondary_y": True}, {"secondary_y": True}],
                [{"secondary_y": False}, {"secondary_y": False}],
                [{"secondary_y": True}, {"secondary_y": False}],
            ],
            subplot_titles=[
                "Electrical response",
                "Temperature-compensated current and resistance",
                "Hydraulic flows",
                "Hydraulic pressures",
                "Temperature and fouling trend",
                "Compensated conductivity channels",
            ],
            vertical_spacing=0.09,
            horizontal_spacing=0.08,
        )

        fig.add_trace(
            go.Scatter(
                x=run_df["Time per step (min)"],
                y=run_df["Ustack"],
                mode="lines",
                name="Ustack [V]",
                line={"color": "#2563eb", "width": 2.5},
            ),
            row=1,
            col=1,
            secondary_y=False,
        )
        fig.add_trace(
            go.Scatter(
                x=run_df["Time per step (min)"],
                y=run_df["Rlinear"],
                mode="lines",
                name="Rlinear",
                line={"color": "#dc2626", "width": 2.5},
            ),
            row=1,
            col=1,
            secondary_y=True,
        )

        fig.add_trace(
            go.Scatter(
                x=run_df["Time per step (min)"],
                y=run_df[CURRENT_25_COL],
                mode="lines",
                name="I25 [A]",
                line={"color": "#059669", "width": 2.5},
            ),
            row=1,
            col=2,
            secondary_y=False,
        )
        fig.add_trace(
            go.Scatter(
                x=run_df["Time per step (min)"],
                y=run_df[RESISTANCE_25_COL],
                mode="lines",
                name="R25 [Ohm]",
                line={"color": "#d97706", "width": 2.5},
            ),
            row=1,
            col=2,
            secondary_y=True,
        )

        for column, label, color in [
            ("FIT2", "FIT2 [l/h]", "#0f766e"),
            ("FIT1", "FIT1 [l/h]", "#7c3aed"),
        ]:
            fig.add_trace(
                go.Scatter(
                    x=run_df["Time per step (min)"],
                    y=run_df[column],
                    mode="lines",
                    name=label,
                    line={"color": color, "width": 2.2},
                ),
                row=2,
                col=1,
            )

        for column, label, color in [
            ("PIT4", "PIT4 [bar]", "#1d4ed8"),
            ("PIT3", "PIT3 [bar]", "#b91c1c"),
        ]:
            fig.add_trace(
                go.Scatter(
                    x=run_df["Time per step (min)"],
                    y=run_df[column],
                    mode="lines",
                    name=label,
                    line={"color": color, "width": 2.2},
                ),
                row=2,
                col=2,
            )

        fig.add_trace(
            go.Scatter(
                x=run_df["Time per step (min)"],
                y=run_df["TIT1"],
                mode="lines",
                name="TIT1 [C]",
                line={"color": "#0891b2", "width": 2.5},
            ),
            row=3,
            col=1,
            secondary_y=False,
        )
        fig.add_trace(
            go.Scatter(
                x=run_df["Time per step (min)"],
                y=run_df["delta_rlinear"],
                mode="lines",
                name="Delta Rlinear",
                line={"color": "#7c2d12", "width": 2.5},
            ),
            row=3,
            col=1,
            secondary_y=True,
        )

        for conductivity_col, color in zip(PLOTLY_CONDUCTIVITY_COLS, conductivity_colors):
            fig.add_trace(
                go.Scatter(
                    x=run_df["Time per step (min)"],
                    y=run_df[conductivity_col],
                    mode="lines",
                    name=conductivity_col,
                    line={"color": color, "width": 2.0},
                ),
                row=3,
                col=2,
            )

        for row in range(1, 4):
            for col in range(1, 3):
                fig.update_xaxes(
                    title_text="Time per step [min]",
                    showgrid=True,
                    gridcolor="#d9e2ec",
                    zeroline=False,
                    row=row,
                    col=col,
                )

        fig.update_yaxes(title_text="Ustack [V]", row=1, col=1, secondary_y=False)
        fig.update_yaxes(title_text="Rlinear", row=1, col=1, secondary_y=True)
        fig.update_yaxes(title_text="I25 [A]", row=1, col=2, secondary_y=False)
        fig.update_yaxes(title_text="R25 [Ohm]", row=1, col=2, secondary_y=True)
        fig.update_yaxes(title_text="Flow [l/h]", row=2, col=1)
        fig.update_yaxes(title_text="Pressure [bar]", row=2, col=2)
        fig.update_yaxes(title_text="TIT1 [C]", row=3, col=1, secondary_y=False)
        fig.update_yaxes(title_text="Delta Rlinear", row=3, col=1, secondary_y=True)
        fig.update_yaxes(title_text="Conductivity (25 C compensated)", row=3, col=2)

        fig.update_layout(
            title={"text": title, "x": 0.5},
            template="plotly_white",
            hovermode="x unified",
            height=1180,
            width=1500,
            legend={
                "orientation": "h",
                "yanchor": "bottom",
                "y": 1.02,
                "xanchor": "left",
                "x": 0.0,
            },
            margin={"t": 130, "l": 70, "r": 70, "b": 70},
        )
        fig.write_html(
            output_dir / filename,
            include_plotlyjs="directory",
            full_html=True,
            config={"responsive": True, "displaylogo": False},
        )
        page_rows.append(
            {
                "run_id": str(run_id),
                "v": f"{run_meta['v']:.1f}",
                "Ce": f"{run_meta['Ce']:.1f}",
                "I": f"{run_meta['I']:.2f}",
                "duration_min": f"{duration_min:.2f}",
                "start_time": _format_timestamp(run_df["date/time"].min()),
                "end_time": _format_timestamp(run_df["date/time"].max()),
                "filename": filename,
            }
        )

    _build_plotly_index(output_dir, page_rows)


def main() -> None:
    """Run the full static analysis pipeline and regenerate figures and tables."""
    sns.set_theme(style="whitegrid", context="talk")
    OUTPUT_DIR.mkdir(exist_ok=True)
    PLOTLY_OUTPUT_DIR.mkdir(exist_ok=True)
    cleanup_outputs(OUTPUT_DIR)

    df = load_data(DATA_PATH)
    df = add_run_id(df)
    all_run_summary = build_all_run_summary(df)
    run_window = prepare_run_window(df, RUN_SELECTIONS)
    full_run_window = prepare_run_window(df, FULL_RUN_SELECTIONS)
    active = prepare_active_runs(df)
    active = add_relative_metrics(active)
    plotly_runs = prepare_plotly_runs(df, PLOTLY_RUN_IDS)
    summary = build_run_summary(active)
    window_summary = build_window_summary(active)
    inverse_current = build_inverse_current_summary(active, summary, window_summary)

    save_profiles(active)
    save_run_summary(summary)
    save_run_duration_status_summary(all_run_summary)
    save_inverse_current_plots(inverse_current)
    save_resistance_25_inverse_current(active)
    save_run_window_overview(run_window, RUN_SELECTIONS)
    save_run_window_overview(full_run_window, FULL_RUN_SELECTIONS)
    save_run_window_hydraulics_overview(run_window, RUN_SELECTIONS)
    save_run_window_hydraulics_overview(full_run_window, FULL_RUN_SELECTIONS)
    export_tables(summary, window_summary, inverse_current)
    save_plotly_run_pages(plotly_runs, PLOTLY_OUTPUT_DIR)

    print("Saved figures:")
    for path in sorted(OUTPUT_DIR.glob("*.png")):
        print(path)
    print("Saved tables:")
    for path in sorted(OUTPUT_DIR.glob("*.csv")):
        print(path)
    print("Saved Plotly pages:")
    for path in sorted(PLOTLY_OUTPUT_DIR.glob("*.html")):
        print(path)


if __name__ == "__main__":
    main()
