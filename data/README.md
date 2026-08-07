# Processed field data

`vof_processed_data.csv` is the processed call-sequence dataset used by
`exploreFieldData.m`. It is the starting point for all empirical analyses
in this repository.

The table contains 6,360 call-to-call observations from 496 sequence
identifiers. The analysis script applies robust sequence-wise screening and
a velocity threshold before fitting or comparing field and simulated
timing.

## Columns

| Column | Description | Unit |
|---|---|---|
| `SeqNum` | Processed sequence identifier | — |
| `Rate` | Call rate for the call-to-call interval | Hz |
| `Velocity` | Estimated bat velocity across the interval | m s⁻¹ |
| `IPI` | Inter-pulse interval | s |
| `Distance` | Estimated displacement across the interval | m |
| `DurationDiff` | Change in call duration between successive calls | s |
| `X`, `Y`, `Z` | Localised call position | m |
| `Timestamp` | Call timestamp within the selected sequence | s |
| `Duration` | Call duration | s |
| `Level` | RMS-derived call-level measure produced from the upstream analyser’s reference channel | dB |

## Use in the repository

`exploreFieldData.m` reads this file with a path resolved from the repository
root. It:

- checks values and flags sequence-wise outliers;
- retains the cleaned observations used for modelling;
- derives timing quantities required by the responsivity comparison;
- generates matched simulations;
- compares displacement, velocity, IPI, call duration, and timing residuals;
- optionally writes `results/field_data_exploration/vof_processed_data_cleaned.csv`.

The source multichannel WAV recordings are not included. The interactive
upstream tools are retained in
[`../field_recording_processing/`](../field_recording_processing/), but the
manual selection and curation leading to this CSV cannot be reproduced as
a fully automated transformation from the material currently available.
