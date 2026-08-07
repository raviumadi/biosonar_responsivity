# Raw field-recording processing

This folder preserves the interactive MATLAB tools used upstream of the
processed field table. They operate on multichannel WAV recordings,
identify calls, allow manual validation, estimate call duration and level,
localise the bat from inter-microphone delays, and derive timestamps,
velocity, and call rate.

Raw recordings are not bundled in this repository. These tools are
therefore optional and are not called by `reproduce_all`.

## Portable entry point

From the repository root:

```matlab
addpath("field_recording_processing")

processedFiles = processFieldRecordings("/path/to/raw/wav", ...
    Mode="regular", ...
    MaxFiles=1);
```

Available modes are:

| Mode | Analyser | Intended use |
|---|---|---|
| `"regular"` | `arrayDataAnalyzer` | Threshold-assisted detection and validation of ordinary approach sequences |
| `"buzzManual"` | `arrayDataAnalyzerwithBuzzManual` | Manual call boundaries for densely packed terminal sequences |
| `"buzzAutomatic"` | `arrayDataAnalyzerwithBuzz` | Legacy threshold-assisted buzz workflow |

The default output location is
`results/field_recording_processing/`. Extracted calls are placed in a
separate subfolder for each source WAV, and successful files are recorded
in a processing log. Regular analyses write `AnalysisTable.mat`; buzz
analyses write `AnalysisTable_Buzz.mat`.

By default, the launcher processes one file and leaves the source WAV
untouched. Set `MaxFiles` to process more files. `MoveSourceAudio=true`
retains the historical behaviour of moving an analysed recording into its
analysis directory and should be used only when that relocation is
intended.

## Included classes

- `arrayDataAnalyzer.m` — regular sequence analysis.
- `arrayDataAnalyzerwithBuzzManual.m` — manual terminal-call boundaries.
- `arrayDataAnalyzerwithBuzz.m` — legacy automatic/threshold-assisted buzz
  variant.
- `AudioFileIterator.m` — legacy datastore iterator retained for provenance;
  the portable launcher uses its own completion log and does not require
  this class.

## Requirements and assumptions

The analysers use Signal Processing Toolbox and Curve Fitting Toolbox
functions. `AudioFileIterator` additionally requires Audio Toolbox.

The microphone geometry and speed of sound used for localisation are
defined in each analyser constructor. These values describe the original
field-array implementation and should be checked before applying the tools
to a different recording system.

The workflow is interactive and requires a MATLAB desktop session for
segment selection, threshold placement, call-boundary selection, and call
validation. It does not provide a fully automated recreation of
`data/vof_processed_data.csv`, which also reflects later manual curation.
