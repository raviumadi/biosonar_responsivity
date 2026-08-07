# Validation

This folder contains numerical and integration checks for the consolidated
responsivity simulator.

## `test_simulateResponsivityCore.m`

Runs smoke tests covering:

- a default first-arrival simulation;
- a multi-target, anchor-switching SSG-like condition;
- optional wingbeat-state metadata;
- call/echo waveform synthesis.

Run it from the repository root:

```matlab
run(fullfile("validation", "test_simulateResponsivityCore.m"))
```

The script raises an assertion if a core structural or numerical sanity
check fails.

## `crossCheckCoreCrVr.m`

Checks the internal agreement between the simulator and the motion-aware
responsivity identities. Three scenarios are supported:

- `"exact"` — deterministic first-arrival case with strict numerical
  assertions;
- `"jitteredBat"` — variable bat velocity;
- `"movingTarget"` — stochastic target motion.

`reproduce_all` selects and runs all three scenarios automatically. When
figure saving is enabled, the outputs are written to
`results/validation_figures/` as PDF and PNG files.
