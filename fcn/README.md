# Helper functions

This folder contains shared utilities used by the simulator, analytical
scripts, and figure-export workflow.

| File | Role | Used by |
|---|---|---|
| `formatLatex.m` | Applies print-size figure dimensions, consistent LaTeX typography, and manuscript-scale font sizes | Figure-generating `explore*.m` scripts, SSG pipeline, and validation figure |
| `exportPaperFigure.m` | Exports the final MATLAB figure as vector PDF and 300-dpi PNG | All saved manuscript figures |
| `reserveLegendMetricMargin.m` | Adds controlled legend padding to compensate for differences between screen and vector-PDF text metrics | Composite plots with dense LaTeX legends, particularly field and SSG analyses |
| `generateVirtualBatCall.m` | Generates an FM sweep with windowing and optional spectral shaping | `synthesiseResponsivityAudio.m` |
| `generateCFBatCall.m` | Generates a CF call with optional Doppler resampling | `synthesiseResponsivityAudio.m` |
| `responsivityRunConfig.m` | Holds the temporary figure/table override used by the master runner | `reproduce_all.m` and the script-based workflows |

## Figure workflow

Plotting scripts normally follow this order:

1. construct all axes, annotations, legends, and colour bars;
2. call `formatLatex` with the intended printed-size preset;
3. apply `reserveLegendMetricMargin` where a boxed legend needs extra vector-export allowance;
4. call `exportPaperFigure`.

Formatting at the intended printed dimensions prevents labels from becoming
disproportionate when a figure is later inserted into a single- or
double-column manuscript layout.

## Path use

Scripts add this folder relative to their own location. No manual `addpath`
configuration is needed when a top-level analysis or `reproduce_all` is
used.
