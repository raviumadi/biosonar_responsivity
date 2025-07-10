# Biosonar Responsivity - Framework with Simulated Bat *Vespertillio numeralis*

This repository contains the MATLAB code and data used in the analyses and simulations for the paper:

**[Biosonar Responsivity Sets the Stage for the Terminal Buzz](https://www.biorxiv.org/content/10.1101/2025.06.16.659925v1)**  

The codebase implements the framework for computing *responsivity*, a new metric that captures the fine-grained temporal precision in echolocation behaviour, and validates it against both simulated and high-resolution field datasets. This toolkit provides scripts to simulate target pursuit sequences, extract call dynamics from real bat recordings, and compute critical behavioural transitions such as the **buzz readiness point** and **reaction time limit** ($$T_{b^*}$$).

## Related Codebase

See [Analyse Responsivity](https://github.com/raviumadi/analyse_responsivity) for a related set of tools. 

---

##  Field Data Analysis Kit `/mat`

### `arrayDataAnalyzer.m`  
Core class for parsing multi-channel audio files from the field array recordings. Includes methods for extracting call timing, durations, computing call rate, estimating velocity, and segmenting for responsivity analyses.

### `arrayDataAnalyzerwithBuzz.m`  
An extended version of `arrayDataAnalyzer` that includes semi-automated buzz detection and validation routines for aligning and extracting detected calls.

### `AudioFileIterator.m`

Utility class that iterates through a set of multi-channel .wav files and applies batch analysis. Used during automated dataset construction.

### `analyzer_files.m`  
Helper script that loads field datasets, applies relevant functions in `arrayDataAnalyzer`, and saves the intermediate output to `.mat` files for further analysis. Run through files while processing data. 

## Data

The outputs of the field data processing are stored in `results` and the processed data are stored in `data`. 

---

## Simulation Scripts `/mat`

### `simulating_sequences.m`

Script to generate synthetic foraging trajectories with randomised prey motion. The resulting call sequences model the bat’s response to variable relative velocity, allowing comparison with field data. Here, used to generate plots. Adopt as fits.

### `foraging_sequence.m`  
Function used within the simulation pipeline to generate an individual prey–bat trajectory, compute call emission timing, and apply the responsivity framework.

---

## Responsivity and Modelling

### `responsivity_curves.m`  
Extracts and plots responsivity metrics across both real and simulated datasets. Outputs aligned timelines with `Tb*`, `Rmax`, and buzz onset.

### `call_duration_model.m`  
Calculates and plots the onset of call duration contraction across sequences. Used to validate that shortening occurs in accordance with decreasing echo delay and increasing call rate.

### `reaction_times.m`

Spatiotemporal parameter plots.

### `delta_plots.m`

Models and call rate against $$k_r$$ and velocity.

### `time_to_intercept.m`

Models the time to intercept based on the relative velocity between the chaser and the target.

## 🧪 Getting Started

Ensure you have MATLAB R2020b or later. Run the following script to add all folders to the path:

```
init
```

---

## Citation

If you use this code in your research, please cite the paper and link to this repository.

---

## License

This work is distributed under the **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)** License.


## Other Publications

Check my [ResearchGate](https://www.researchgate.net/profile/Ravi-Umadi-3) or [ORCID](https://orcid.org/0000-0003-3867-1769) for a full list of my research work. 

## Contact

Drop by my personal website [biosonix.io](https://biosonix.io) and drop a message if you would like to collaborate or need assistance with the code and development. 

---