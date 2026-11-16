# MIAD paper code and data

This folder contains the MATLAB workflows and experimental materials supporting
the MIAD manuscript.

## Repository layout

- `code/run_simulation.m` reproduces the Monte Carlo comparison and exports
  the RMSE figure and numerical results.
- `code/run_experiment.m` processes the three stationary intervals in the
  4 September 2025 IMU recording and exports the attitude comparison.
- Other MATLAB files in `code/` implement the MIAD candidates, candidate
  selection, uncertainty analysis, attitude conversion, and comparison
  algorithms.
- `code/LICENSE` contains the MIT license for project-authored MATLAB files
  marked with `SPDX-License-Identifier: MIT`.
- `data/experimental/2025-09-04/raw/` contains the unmodified sensor logs.
- `data/experimental/2025-09-04/photos/` contains photographs documenting
  the three test positions.
- `results/` is populated by the run scripts and is not required as input.

## Notation

The code uses the manuscript's vector notation: `a` and `b` are the observed
sensor-frame vectors; `p` and `q` are the corresponding reference-frame unit
vectors; `i` is the index of the unobserved component of `a`; and `j` is the
index of the selected primary component of `b`. In the scripts these appear
as `reference_p`, `reference_q`, `unknown_index_i`, and
`observed_index_j`. The candidate arrays are named `a_candidates` and
`b_candidates`; `sigma_a_i1`, `sigma_a_i2`, and `sigma_b_j` denote the
standard deviations of the three primary measurements.

## Requirements and running

Tested with MATLAB R2024b. The workflows use base MATLAB functionality and
write their outputs beneath `results/`.

From MATLAB, run either script by full path, for example:

```matlab
run(fullfile(pwd,'code','run_simulation.m'))
run(fullfile(pwd,'code','run_experiment.m'))
```

Run the commands from this folder. The scripts set their own code and data
paths. The simulation generates 50 test cases with 1,000 noise realizations per case.
The experimental script uses the recorded time windows 20–60 s,
120–160 s, and 240–280 s, measured from the first timestamp.

## Experimental data

The ZIMU070 text log contains timestamp in milliseconds, three gyroscope
channels in degrees per second, and three accelerometer channels in metres
per second squared. The ER-3MG-01 text log contains timestamp in milliseconds
and three gyroscope channels in degrees per second. The experiment script
applies the channel mapping and unit conversions used for the analysis. It
then applies a right-handed 90-degree in-plane body-axis rotation so that the
selected primary gyroscope measurement is indexed as `b_1`, matching
`(i,j)=(3,1)` in the manuscript. This maps the second recorded ZIMU070 gyro
channel to `b_1`; verify the physical channel labels against the manuscript's
sensor-axis description before submission, since the raw log has no header.
The recorded positions were reached by successive 90-degree clockwise
rotations from Position 1 to Position 3, followed by a counter-clockwise
return to Position 2. Recording stopped during the final Position 2 interval
after the FOG acquisition interface stopped recording; no samples for that
final interval are used in the analysis.

The experimental estimates should be checked against the manuscript's
reported table after each change to the axis mapping, time windows, or sensor
calibration.

## Generated files

The simulation writes `figure_1_attitude_rmse.pdf`,
`simulation_results.csv`, `simulation_diagnostics.csv`, and
`simulation_results.mat` to `results/simulation/`. It also exports
`norm_violation_results.csv`, `negative_delta_results.csv`,
`category1_misdetection_results.csv`, and
`category2_misdetection_results.csv`, containing the data used in the simulation. The diagnostics CSV contains the candidate-selection
wrong rates and mean error values used for the degeneracy analysis. The
experiment writes `experiment_results.csv` and `experiment_results.mat` to
`results/experiment/`.
