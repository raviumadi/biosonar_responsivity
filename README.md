# The Temporal Logic of Bat Echolocation: A Responsivity Theory of Call Timing Patterns

This repository contains the MATLAB implementation, simulations, validation
checks, processed field data, and figure-generating analyses associated with
the manuscript **“The Temporal Logic of Bat Echolocation: A Responsivity
Theory of Call Timing Patterns”** by Ravi Umadi and Uwe Firzlaff.

The manuscript was previously titled
[*Biosonar Responsivity Sets the Stage for the Terminal Buzz*](https://www.biorxiv.org/content/10.1101/2025.06.16.659925v2).
The linked bioRxiv version preserves the earlier stage of the work and its
original title; the current title reflects the expanded theory and analyses.

The responsivity framework describes how the timing of one active-sensing
event can be decomposed into acoustic information acquisition and a
behavioural response opportunity. The code develops this event-level
description into call-by-call simulations and asks how call timing changes
with acoustic delay, relative motion, the effective acoustic anchor,
responsivity coefficient, call-duration constraints, and stochastic
variation. Downstream analyses examine terminal timing, buzz readiness,
sonar strobe groups (SSGs), and compatibility with field recordings from
free-flying bats.

## Quick start

Clone the repository, open MATLAB in its root directory, and run:

```matlab
summary = reproduce_all;
```

This regenerates the analytical figures and tables under `results/`, runs
the SSG simulation, and executes all three core-validation scenarios. The
returned table reports whether each component completed and how long it
took.

For a test run that displays analyses without rewriting outputs:

```matlab
summary = reproduce_all( ...
    SaveFigures=false, ...
    SaveTables=false, ...
    RunSSG=false);
```

The SSG parameter sweep is the longest component. It can be omitted during
development with `RunSSG=false`. Validation can likewise be controlled with
`RunValidation=false`; `StopOnError=false` allows the remaining components
to continue if one analysis fails. By default, `reproduce_all` keeps every
generated figure open for inspection. Use `KeepFigures=false` to close the
previous analysis figure before each new component and reduce desktop
clutter.

## Requirements

The package has been tested with **MATLAB R2025b**.

- MATLAB
- Statistics and Machine Learning Toolbox for the statistical analyses
- Signal Processing Toolbox for call/echo synthesis and its smoke test

The optional raw-recording workflow also uses Signal Processing Toolbox and
Curve Fitting Toolbox functions. The legacy `AudioFileIterator` utility
uses Audio Toolbox, but it is not required by the portable raw-processing
entry point or by `reproduce_all`.

## Repository structure

```text
.
├── reproduce_all.m                 master reproduction entry point
├── simulateResponsivityCore.m      unified call-by-call simulator
├── synthesiseResponsivityAudio.m   optional call and echo synthesis
├── markSSGPatterns.m               downstream SSG detector
├── explore*.m                      manuscript analyses and figures
├── runCoreSSGAnalysisPipeline.m    multi-condition SSG analysis
├── data/                           processed field dataset
├── fcn/                            plotting and audio helper functions
├── field_recording_processing/     optional upstream WAV analysis
├── validation/                     numerical checks and smoke tests
└── results/                        figures and derived tables
```

Each principal subfolder contains a short README describing its contents:

- [`data/README.md`](data/README.md)
- [`fcn/README.md`](fcn/README.md)
- [`field_recording_processing/README.md`](field_recording_processing/README.md)
- [`validation/README.md`](validation/README.md)
- [`results/README.md`](results/README.md)

## Main analytical workflows

| Script | Purpose | Principal outputs |
|---|---|---|
| `exploreCallRateDistnaceDynamics.m` | Analytical call-rate–distance curves across responsivity coefficients | `explore_Cr_vs_distance_by_kr` |
| `exploreCallDurationByKr.m` | Call-duration contraction under the no-overlap constraint | `explore_Tc_vs_distance_by_kr` |
| `exploreIPIDistancePhiMotile.m` | Effects of target motion, bat-velocity variation, and information fraction on IPI, call rate, and duration | `explore_IPI_vs_distance_phi_motile` and summary CSV files |
| `exploreTaDistanceDynamics.m` | Acquisition and behavioural timing across responsivity, information, and velocity sweeps | `explore_Ta_vs_distance_dynamic` and summary CSV files |
| `exploreResponsivityCurve.m` | Responsivity curves, buzz readiness, and terminal timing proxies | Buzz-readiness figures and sequence-level summaries |
| `exploreTerminalTbWindow.m` | Minimum behavioural interval and the duration/call count of the terminal spatial window | `explore_shortest_Tb_at_dstop` |
| `exploreFieldData.m` | Quality control, velocity-matched field–simulation comparisons, seeded duration-envelope analysis, and conditional acquisition-window diagnostics | Field–simulation composite figures, profile summaries, and reproducibility caches |
| `runCoreSSGAnalysisPipeline.m` | SSG emergence across single/multiple and stationary/moving target conditions | SSG event tables and composite figure |

The scripts may also be run individually. Their local `savePlots`,
`saveStats`, or `writeTables` switches are intended for figure development.
When called through `reproduce_all`, the master options override these local
switches without changing the scripts. Individual scripts close stale
figures when run directly; the master runner suppresses that behaviour while
`KeepFigures=true`.

## Using the core simulator directly

The simulator accepts parameter and option structures and returns a
call-level table together with the resolved configuration and summary:

```matlab
params = struct( ...
    'kr', 5, ...
    'initialDistance_m', 5, ...
    'initialBatSpeed_m_s', 5);

opts = struct( ...
    'geometryMode', "3D", ...
    'timingMode', "motionAware", ...
    'rngSeed', 1);

result = simulateResponsivityCore(params, opts);
head(result.calls)
```

`simulateResponsivityCore` generates event timing, trajectories, anchor
state, call-duration metadata, and optional wingbeat state. It deliberately
does not label buzzes, SSGs, synchrony, or behavioural feasibility; those
interpretations belong to downstream analyses.

Optional acoustic demonstrations can be generated from a simulation:

```matlab
audio = synthesiseResponsivityAudio(result);
soundsc(audio.stereo, audio.fs);
```

## Field data and upstream recording analysis

The manuscript field analysis begins with
[`data/vof_processed_data.csv`](data/vof_processed_data.csv). This processed
table is included so that the field figures and statistics can be
regenerated without access to the original multichannel recordings.

The optional interactive workflow that detects, validates, localises, and
measures calls from raw WAV files is preserved in
[`field_recording_processing/`](field_recording_processing/). Raw
recordings are not bundled here, and this upstream workflow is therefore
not part of `reproduce_all`. The earlier public repository contains a
subset of example recordings in its
[pre-consolidation revision](https://github.com/raviumadi/biosonar_responsivity/tree/135dfea).

The exact curation that produced the supplied processed CSV included manual
sequence selection and is not represented as a fully automated
raw-to-CSV transformation. The included WAV tools document and preserve
the upstream call extraction and localisation procedure; the supplied CSV
is the reproducible starting point for the analyses reported here.

## Outputs and reproducibility

All load and save paths used by the manuscript workflows are resolved from
the repository location. The code does not depend on the caller’s current
working directory.

- Figures are exported as vector PDF and 300-dpi PNG files.
- Derived statistics and call/event tables are written as CSV files.
- Stochastic manuscript analyses use fixed random seeds.
- The field comparison caches compatible simulation tables and automatically
  rebuilds them when its analysis parameters change.
- Existing outputs with the same names are overwritten when saving is
  enabled.
- The raw-recording launcher does not move source WAV files unless
  `MoveSourceAudio=true` is explicitly supplied.

See [`results/README.md`](results/README.md) for the output-to-script map.

## Validation

Run the standalone smoke tests with:

```matlab
run(fullfile("validation", "test_simulateResponsivityCore.m"))
```

The numerical validation compares stored simulator quantities with the
motion-aware responsivity identities in deterministic, jittered-bat, and
moving-target conditions. `reproduce_all` runs all three conditions.

## Citation

If you use this framework or code, please cite the latest version of the
associated manuscript:

> Umadi, R. & Firzlaff, U. *The Temporal Logic of Bat Echolocation: A
> Responsivity Theory of Call Timing Patterns*. bioRxiv (2026). [doi: 10.1101/2025.06.16.659925](https://doi.org/10.1101/2025.06.16.659925)

The preceding bioRxiv version remains available under the former title for
version-level traceability:

> Umadi, R. & Firzlaff, U. *Biosonar Responsivity Sets the Stage for the
> Terminal Buzz*. bioRxiv, version 2 (2025).
> <https://www.biorxiv.org/content/10.1101/2025.06.16.659925v2>

The citation metadata will be updated when the revised bioRxiv version or
the final article receives its definitive publication details.

## Licence

This repository is distributed under the
[Creative Commons Attribution-NonCommercial 4.0 International licence](LICENSE.md).
