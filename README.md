# Electrodialysis Data Processing And Visualization

This folder contains a static analysis workflow for the electrodialysis fouling dataset stored in [dataRaw.csv](dataRaw.csv).

The current focus is on active runs at:

- `2.5 cm/s`
- `3.0 cm/s`
- `3.5 cm/s`

using Seaborn, with layouts chosen to support engineering comparison of current, voltage, resistance, flow, pressure, and temperature behavior.

## Main Files

- [visualization.py](visualization.py): loads the CSV, builds derived metrics, and writes figures and tables
- [dataRaw.csv](dataRaw.csv): optional local working copy of the raw time-series dataset
- [requirements.txt](requirements.txt): Python dependencies for the workflow

## External Sources

- Article page: <https://www.sciencedirect.com/science/article/pii/S2352340920306570>
- Dataset record: <https://zenodo.org/records/3551928>
- Dataset DOI: <https://doi.org/10.5281/zenodo.3551928>

For a public GitHub repository, these external links are the preferred references.
If you do not want to commit local source copies, keep the links above and place the files locally only when running the analysis.

## License Notes

This repository can contain materials under different licenses.

- Code written for this analysis can be licensed separately by the repository maintainer, for example under `MIT`
- The source dataset file [dataRaw.csv](dataRaw.csv), if included in the repository, should remain attributed to its original source and treated as `Creative Commons Attribution 4.0 International (CC BY 4.0)`
- Top-level license files: [LICENSE](LICENSE) for repository code and [DATA_LICENSE.md](DATA_LICENSE.md) for the dataset notice

If `dataRaw.csv` is published in the repo, keep the dataset attribution with:

- online data source: <https://zenodo.org/records/3551928>
- Zenodo record: <https://zenodo.org/records/3551928>
- DOI: <https://doi.org/10.5281/zenodo.3551928>
- license: `CC BY 4.0`

Ready-to-use attribution sentence:

`Dataset source: Zenodo, https://zenodo.org/records/3551928, DOI: 10.5281/zenodo.3551928, licensed under CC BY 4.0.`

Do not describe the whole repository as `MIT only` when the CSV dataset is included.

## Documentation

- [DATA_MANIPULATION_WORKFLOW.md](DATA_MANIPULATION_WORKFLOW.md): pipeline explanation and formulas
- [DELTA_RLINEAR_EXPLANATION.md](DELTA_RLINEAR_EXPLANATION.md): explanation of `delta_rlinear`

## Current Outputs

Generated figures and tables are written to [outputs](outputs).
Interactive Plotly pages for selected runs are written to [outputs_plotly](outputs_plotly).

Key figures:

- [voltage_article_style.png](outputs/voltage_article_style.png)
- [resistance_25_article_style.png](outputs/resistance_25_article_style.png)
- [resistance_25_vs_inverse_current.png](outputs/resistance_25_vs_inverse_current.png)
- [rlinear_vs_inverse_current.png](outputs/rlinear_vs_inverse_current.png)
- [delta_rlinear_profiles.png](outputs/delta_rlinear_profiles.png)
- [norm_delta_rlinear_profiles.png](outputs/norm_delta_rlinear_profiles.png)
- [run_duration_status_scatter.png](outputs/run_duration_status_scatter.png)
- [run_1_14_29_45_current_resistance_overview.png](outputs/run_1_14_29_45_current_resistance_overview.png)
- [run_1_14_29_45_diluate_concentrate_overview.png](outputs/run_1_14_29_45_diluate_concentrate_overview.png)
- [run_1_45_current_resistance_overview.png](outputs/run_1_45_current_resistance_overview.png)
- [run_1_45_diluate_concentrate_overview.png](outputs/run_1_45_diluate_concentrate_overview.png)

Key exported tables:

- [run_summary.csv](outputs/run_summary.csv)
- [window_summary.csv](outputs/window_summary.csv)
- [inverse_current_summary.csv](outputs/inverse_current_summary.csv)

Plotly exports:

- [index.html](outputs_plotly/index.html)
- `run_15_interactive.html`
- `run_17_interactive.html`
- `run_19_interactive.html`
- `run_21_interactive.html`

## How To Run

Install dependencies:

```powershell
python -m pip install -r requirements.txt
```

If `dataRaw.csv` is not present locally, download it from the Zenodo dataset record and place it in this directory.

Run the workflow from this directory:

```powershell
python visualization.py
```

## Current Scope

The present workflow:

- keeps only active runs with positive current for the main comparison plots
- uses `25 deg C` normalized current with `1 % / deg C` compensation

$$
I_{25} = \frac{I}{1 + 0.01(T-25)}
$$

- uses `25 deg C` normalized apparent resistance

$$
R_{25} = \frac{U_{\mathrm{stack}}}{I_{25}}
$$

- computes conductivity-normalized channels with `2 % / deg C` compensation
- keeps voltage plots in measured form
- keeps `Rlinear` as a separate processed resistance-like trend signal
- includes stitched cumulative-time overview plots for selected run blocks and for the full `1-45` run range
- exports standalone Plotly HTML pages for runs `15`, `17`, `19`, and `21`

It does not attempt a full mechanistic limiting-current model or an interactive dashboard.
