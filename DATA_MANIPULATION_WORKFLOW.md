# Data Manipulation Workflow

This note explains how the raw file [dataRaw.csv](c:/Users/tomas/OneDrive/python/robodreams/rd-python-industrial-automation-v2/tests/leson3/ukol3/dataRaw.csv) is transformed inside [visualization.py](c:/Users/tomas/OneDrive/python/robodreams/rd-python-industrial-automation-v2/tests/leson3/ukol3/visualization.py) before plotting.

The goal is simple:

1. isolate comparable electrodialysis runs
2. compare `2.5`, `3.0`, and `3.5 cm/s`
3. show voltage and resistance trends in a layout similar to the article
4. derive a few resistance-growth indicators that help interpret fouling behavior

> [!NOTE]
> This document describes the code that is currently in the repository.
> It is not a full scientific interpretation of the process.
> It is a transparent record of what was done to the data.

## Input Data

The CSV contains time-series measurements such as:

- `date/time`
- `Time (h)`
- `Time per step (min)`
- `Time per cicle (h)`
- `Resistance (Ohm)`
- `Ustack`
- `Rlinear`
- `Ce`
- `v`
- `I`
- pressure, pH, temperature, conductivity, and flow columns

The core variables used in the current workflow are:

- $v$: linear flow velocity
- $C_e$: inlet salt concentration
- $I$: applied current
- $I_{25}$: current normalized to `25 ^\circ C`
- $U_{\mathrm{stack}}$: stack voltage
- $R_{25}$: apparent resistance normalized to `25 ^\circ C`
- $R_{\mathrm{linear}}$: processed resistance-like signal from the dataset

The current normalization used in the plots is:

$$
I_{25} = \frac{I}{1 + 0.01(T-25)}
$$

The apparent resistance used in the plots is then defined as:

$$
R_{25} = \frac{U_{\mathrm{stack}}}{I_{25}}
$$

and the compensated power is:

$$
P_{25} = U_{\mathrm{stack}} \cdot I_{25}
$$

The conductivity compensation used for the `CO*` channels is:

$$
\kappa_{25} = \frac{\kappa_T}{1 + 0.02(T-25)}
$$

> [!IMPORTANT]
> The CSV also contains a column named `Resistance (Ohm)`, but the earlier sanity check showed that it does not match the direct quotient $U_{\mathrm{stack}} / I$ for the `0.2 M` runs.
> Because of that mismatch, the current workflow uses the compensated apparent resistance $R_{25}$ as the main resistance quantity for interpretation and plotting.

## Step 1: Load And Parse

The CSV is read with `pandas.read_csv`, and the `date/time` column is parsed into timestamps.

No rows are dropped at this stage.

## Step 2: Split The Dataset Into Runs

The file contains multiple experiments placed one after another.
To avoid drawing one long line across unrelated experiments, the script creates a `run_id`.

A new run starts when at least one of the following happens:

- `Time per step (min)` goes backward
- `v` changes
- `Ce` changes
- `I` changes

In compact form, a new run is triggered when

$$
\mathrm{reset}(t) =
\left[\Delta t_{\mathrm{step}} < 0\right]
\lor
\left[\Delta v \neq 0\right]
\lor
\left[\Delta C_e \neq 0\right]
\lor
\left[\Delta I \neq 0\right]
$$

and the run index is the cumulative count of resets:

$$
\mathrm{run\_id}(k) = 1 + \sum_{j=1}^{k} \mathrm{reset}(j)
$$

> [!IMPORTANT]
> This is an operational segmentation rule.
> It is meant to identify contiguous experimental blocks in the CSV, not to re-derive the original lab protocol from first principles.

## Step 3: Keep Only The Relevant Operating Window

The current plotting workflow keeps only:

- active electrodialysis periods: $I > 0$
- requested flow velocities: $v \in \{2.5,\ 3.0,\ 3.5\}$

So the filtered dataset is:

$$
\mathcal{D}_{\mathrm{active}} =
\{ x \in \mathcal{D} \mid I(x) > 0,\; v(x) \in \{2.5, 3.0, 3.5\} \}
$$

This removes the polarity-reversal cleaning segments with negative current from the default comparison plots.

> [!TIP]
> This is why the main figures focus on fouling buildup during the active cycle rather than on the cleaning cycle.

## Step 4: Create Human-Readable Labels

For faceting and legends, the script creates label columns:

- `v_label`, for example `3.0 cm/s`
- `Ce_label`, for example `0.1 M`
- `I_label`, for example `1.52 A`

This does not change the data numerically.
It only makes the plots readable.

## Step 5: Derive Relative Resistance Metrics

For each run, the first active sample is used as the baseline.

If the baseline value is at time $t_0$, then:

$$
R_{\mathrm{linear},0} = R_{\mathrm{linear}}(t_0)
$$

and the run-relative change is:

$$
\Delta R_{\mathrm{linear}}(t) =
R_{\mathrm{linear}}(t) - R_{\mathrm{linear}}(t_0)
$$

The normalized version of `delta_rlinear` is:

$$
\Delta R_{\mathrm{linear,norm}}(t) =
\frac{R_{\mathrm{linear}}(t) - R_{\mathrm{linear}}(t_0)}
{R_{\mathrm{linear}}(t_0)}
$$

These transformations make different runs comparable even if they start from different baseline resistance levels.

## Step 6: Build Run-Level Summaries

Each run is summarized by taking:

- first `v`, `Ce`, `I`
- first `I25`
- run duration
- initial and final `Rlinear`
- initial and final `R25`
- final `Ustack`

The main run-level fouling indicator is:

$$
\Delta R_{\mathrm{linear,run}} =
R_{\mathrm{linear}}(t_{\mathrm{end}}) - R_{\mathrm{linear}}(t_0)
$$

This is exported to [run_summary.csv](c:/Users/tomas/OneDrive/python/robodreams/rd-python-industrial-automation-v2/tests/leson3/ukol3/outputs/run_summary.csv).

## Step 7: Build Early/Late Window Summaries

To avoid depending only on one endpoint sample, the script also computes windowed medians for each run:

- early window: first `10` minutes
- late window: last `10` minutes of the run

For example:

$$
U_{\mathrm{early}} = \mathrm{median}\left(U_{\mathrm{stack}}(t) \;|\; t \le 10\ \mathrm{min}\right)
$$

$$ 
U_{\mathrm{late}} = \mathrm{median}\left(U_{\mathrm{stack}}(t) \;|\; t \ge t_{\max} - 10\ \mathrm{min}\right)
$$

and similarly for the compensated apparent resistance:

$$
R_{25,\mathrm{early}} = \mathrm{median}\left(R_{25}(t) \;|\; t \le 10\ \mathrm{min}\right)
$$

$$
R_{25,\mathrm{late}} = \mathrm{median}\left(R_{25}(t) \;|\; t \ge t_{\max} - 10\ \mathrm{min}\right)
$$

This summary is exported to [window_summary.csv](c:/Users/tomas/OneDrive/python/robodreams/rd-python-industrial-automation-v2/tests/leson3/ukol3/outputs/window_summary.csv).

## Step 8: Plotting Strategy

The main static figures use:

- `Time per step (min)` on the x-axis
- one row per flow velocity
- one column per concentration
- line color for current level
- one line per `run_id`

In Seaborn terms, the important choice is:

- `units="run_id"`
- `estimator=None`

That means each run is drawn directly and Seaborn does not average separate runs together.

> [!NOTE]
> This matters because averaging across runs can hide abrupt changes, truncation, or saturation behavior.

## Current Outputs

The most important current outputs are:

- [voltage_article_style.png](c:/Users/tomas/OneDrive/python/robodreams/rd-python-industrial-automation-v2/tests/leson3/ukol3/outputs/voltage_article_style.png)
- [resistance_25_article_style.png](c:/Users/tomas/OneDrive/python/robodreams/rd-python-industrial-automation-v2/tests/leson3/ukol3/outputs/resistance_25_article_style.png)
- [resistance_25_vs_inverse_current.png](c:/Users/tomas/OneDrive/python/robodreams/rd-python-industrial-automation-v2/tests/leson3/ukol3/outputs/resistance_25_vs_inverse_current.png)
- [rlinear_profiles.png](c:/Users/tomas/OneDrive/python/robodreams/rd-python-industrial-automation-v2/tests/leson3/ukol3/outputs/rlinear_profiles.png)
- [rlinear_vs_inverse_current.png](c:/Users/tomas/OneDrive/python/robodreams/rd-python-industrial-automation-v2/tests/leson3/ukol3/outputs/rlinear_vs_inverse_current.png)
- [delta_rlinear_profiles.png](c:/Users/tomas/OneDrive/python/robodreams/rd-python-industrial-automation-v2/tests/leson3/ukol3/outputs/delta_rlinear_profiles.png)
- [norm_delta_rlinear_profiles.png](c:/Users/tomas/OneDrive/python/robodreams/rd-python-industrial-automation-v2/tests/leson3/ukol3/outputs/norm_delta_rlinear_profiles.png)

## Interpretation Boundaries

> [!CAUTION]
> The transformations in this workflow help compare runs, but they do not remove all physical ambiguity.
> A rise in voltage or resistance can reflect a mixture of:
> - fouling
> - concentration polarization
> - hydrodynamic effects
> - temperature behavior
> - approach to limiting current

So the workflow should be read as:

- a transparent data-processing pipeline
- a plotting pipeline for comparison
- a first-pass interpretation aid

and not as a complete mechanistic model.

## Short Version

> [!NOTE]
> The data manipulation pipeline does five main things:
> 1. splits the CSV into separate runs
> 2. keeps only active runs at `2.5`, `3.0`, and `3.5 cm/s`
> 3. computes `25 °C` normalized current, resistance, power, and conductivity using `1 % / °C` for current and `2 % / °C` for conductivity
> 4. builds run summaries and early/late window summaries
> 5. plots voltage and resistance in a consistent article-style layout with 15-minute time ticks
