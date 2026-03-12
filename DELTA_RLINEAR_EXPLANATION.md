# Delta Rlinear Explanation

> [!NOTE]
> `delta_rlinear` is a simple per-run fouling-growth metric.
> I compute it as the increase of `Rlinear` relative to the start of the active electrodialysis run:
>
> $$
> \Delta R_{\mathrm{linear}}(t) = R_{\mathrm{linear}}(t) - R_{\mathrm{linear}}(t_0)
> $$
>
> where $t_0$ is the first sample of the active run (`I > 0`).

## 1. What `Rlinear` means here

> [!NOTE]
> In the CSV there are two resistance-related columns:
> - `Resistance (Ohm)`
> - `Rlinear`
>
> In the current analysis, `Rlinear` is treated as the processed resistance-like signal for trend comparison.
> For the apparent-resistance plots, I now use a `25 °C` normalized quantity based on a compensated current
>
> $$
> I_{25} = \frac{I}{1 + 0.01(T-25)}
> $$
>
> $$
> R_{25} = \frac{U_{\mathrm{stack}}}{I_{25}}
> $$
>
> instead of relying on the CSV `Resistance (Ohm)` column directly.

> [!IMPORTANT]
> This is partly an inference from the dataset structure and the paper context.
> The article explains that the data were processed to reduce temperature effects and then filtered/resampled.
> So `Rlinear` is useful as a cleaner trend variable, while `Resistance (Ohm)` is still worth plotting as the direct measured quantity.

## 2. Why I used a delta instead of the absolute value

If we only look at the absolute resistance level, we mix together two things:

1. the starting condition of the run
2. the resistance growth during the run

For fouling interpretation, the second part is usually what matters more.

So I subtract the starting value of each active run:

$$
\Delta R_{\mathrm{linear}}(t) = R_{\mathrm{linear}}(t) - R_{\mathrm{linear},0}
$$

with

$$
R_{\mathrm{linear},0} = R_{\mathrm{linear}}(t_0)
$$

This forces every run to start at zero and makes the curves comparable as "growth from that run's own baseline".

> [!WARNING]
> If one run starts at `3.8` and ends at `18.6`, and another starts at `4.3` and ends at `5.8`, the absolute end values tell only part of the story.
>
> The deltas make the difference explicit:
>
> $$
> 18.6 - 3.8 \approx 14.8
> $$
>
> $$
> 5.8 - 4.3 \approx 1.5
> $$
>
> So the first run accumulated much more resistance growth during the active period.

## 3. How it is computed in the code

In [visualization.py](c:/Users/tomas/OneDrive/python/robodreams/rd-python-industrial-automation-v2/tests/leson3/ukol3/visualization.py), the logic is:

1. split the dataset into individual runs using resets in time and changes in `v`, `Ce`, or `I`
2. keep only active electrodialysis runs with positive current
3. for each `run_id`, take the first `Rlinear` value as the run baseline
4. subtract that baseline from every later point in the same run

Mathematically:

$$
\Delta R_{\mathrm{linear},i}(t) = R_{\mathrm{linear},i}(t) - R_{\mathrm{linear},i}(t_{0,i})
$$

where $i$ is the run index.

The normalized version is:

$$
\Delta R_{\mathrm{linear},\mathrm{norm}}(t) =
\frac{R_{\mathrm{linear}}(t) - R_{\mathrm{linear}}(t_0)}{R_{\mathrm{linear}}(t_0)}
$$

This tells you the relative increase instead of the absolute increase.

## 4. Physical interpretation

> [!tip]
> A rising `delta_rlinear` means the stack is becoming harder to drive at the same imposed current.
> In practical electrodialysis terms, that usually indicates increasing transport resistance, commonly consistent with fouling, concentration polarization, or both.

If current is fixed and voltage increases over time, then the apparent resistance also increases:

$$
R \sim \frac{U}{I}
$$

So when `delta_rlinear` rises strongly during an active cycle, the system is drifting away from its initial clean or less-fouled state.

> [!warning]
> `delta_rlinear` is a useful operational indicator, but it is not a pure "fouling thickness" measurement.
> It can contain contributions from:
> - membrane fouling
> - concentration polarization
> - hydrodynamic effects
> - temperature-related behavior not fully removed by processing
> - approaching a limiting-current regime

## 5. What we observed in your data

From the current run summaries:

- `3.0 cm/s`, `0.1 M` shows very large growth, roughly $\Delta R_{\mathrm{linear}} \approx 14.4$ to `16.5`
- `3.5 cm/s`, `0.1 M` increases more progressively, from about `1.5` at low current to about `14.1` at the highest current
- `0.2 M` cases are much smaller overall, except for a clear jump at higher current in `3.0 cm/s`

> [!NOTE]
> The main value of `delta_rlinear` is that it separates "how much the run changed" from "where the run started".
> That is why it is useful for comparing fouling trajectories between `2.5`, `3.0`, and `3.5 cm/s`.

## 6. Why I also use calculated resistance

`delta_rlinear` is good for interpretation, but it is also useful to compare it against an apparent resistance built from the electrical definition after temperature normalization of current:

$$
I_{25} = \frac{I}{1 + 0.01(T-25)}
$$

$$
R_{25} = \frac{U_{\mathrm{stack}}}{I_{25}}
$$

This is now the main resistance definition used in the resistance-based plots because it is transparent and internally consistent with the measured voltage, applied current, and the chosen engineering compensation.

So the two views answer different questions:

- `R_{25}`: what the temperature-normalized apparent electrical resistance looks like
- `delta_rlinear`: how much resistance accumulated relative to the start of each run

## 7. Human summary

> [!WARNING]
> Read `delta_rlinear` as:
> "How much extra resistance did this run build up since the active cycle started?"
>
> If the curve stays near zero, the system stayed close to its initial state.
> If the curve climbs fast, the run is accumulating resistance quickly.
> That makes it a practical fouling-growth indicator for comparing operating conditions.
