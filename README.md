# Venus Spin–Climate Evolution Model (Venus-SCEM)

This repository contains the code needed to run the Venus Spin–Climate Evolution Model (Venus-SCEM) and reproduce the principal results and figures from:

> “No Recent Long-Lived Habitable States on Venus Due to Spin–Climate Interaction”

The model code, analysis scripts, and plotting scripts are hosted in this GitHub repository. The full simulation outputs used in the paper (`all_solutions.pkl`) are archived separately on Zenodo because of their large file sizes.

**Simulation output archive:** ZENODO_DOI_OR_URL

## Repository structure

### `model/`

Core model code and input data.

- **`run_Venus_evolution.py`**  
  Main model script. It reads the sensitivity-test parameters specified in `venus_run_input`, runs the simulations, and sorts the solutions according to their final spin states.

- **`venus_run_input`**  
  Input parameters for the main model script. See `run_Venus_evolution.py` for descriptions of the individual parameters.

- **`run_Venus_many.sh`**  
  Shell script used to submit simulations as batch jobs.

- **`concatenate_sorting.py`**  
  Merges the batch outputs into a single `all_solutions.pkl` file and summarizes the final-state distribution. This script is called automatically after all batch jobs finish successfully. After merging, it removes the individual batch directories.

#### `model/input_data/`

- **`albedo_Yang_matrix.npy`**  
  NumPy array containing planetary albedo as a function of stellar flux and rotation rate, based on Yang et al. (2013). The model interpolates these data to estimate cloud-driven albedo in the habitable state.

- **`Q_withIT_330.npy`**  
  Tidal quality factors derived from the ocean-tide results of Green et al. (2019). These values are interpolated in the sensitivity test that includes ocean tidal dissipation.

- **`sig_330_Green.npy`**  
  Corresponding tidal forcing frequencies for the Green et al. (2019) data.

## Simulation outputs

The complete model outputs used to produce the results in the paper are archived on Zenodo:

**ZENODO_DOI_OR_URL**

Each model ensemble produces an `all_solutions.pkl` file containing the merged simulation results. These files are not stored in the GitHub repository because of their large size.

The Zenodo archive contains the `all_solutions.pkl` files for the baseline simulation and all sensitivity tests described below. To reproduce the figures without rerunning the full ensembles, download the archived outputs from Zenodo and place each `all_solutions.pkl` file in its corresponding directory in this repository.

For example:

```text
baseline/
└── all_solutions.pkl

sensitivity_tests/
├── nu_1e-4/
│   └── all_solutions.pkl
├── nu_10/
│   └── all_solutions.pkl
├── Q_100/
│   └── all_solutions.pkl
└── ...
```

## Running the model

Submit an ensemble of simulations using:

```bash
./run_Venus_many.sh <case_name> <number_of_batches> <runs_per_batch>
```

For example:

```bash
./run_Venus_many.sh baseline 50 1000
```

This command creates a case named `baseline` consisting of 50 batches with 1,000 simulations each, for a total of 50,000 simulations.

> **Note:** `run_Venus_many.sh` is currently configured for the Harvard FAS Research Computing cluster. Users will need to modify the batch-submission settings to run the model on another cluster or local machine.

## Baseline results

### `baseline/`

The baseline ensemble consists of 50,000 simulations and is used to produce Figures 2–4 of the main text.

- **`all_solutions.pkl`**  
  Merged output from the baseline ensemble. Because of its size, this file is archived on Zenodo rather than tracked in this GitHub repository. Download it from the Zenodo archive and place it in `baseline/` to reproduce the analysis and figures.

- **`concatenate_all.out`**  
  Text output containing a summary of the final spin-state distribution.

## Main-text figures

### `figures/`

- **`constants_functions.ipynb`**  
  Defines constants and functions used by the plotting notebooks.

- **`Figure_2.ipynb`**  
  Reproduces Figure 2 of the main text using the baseline `all_solutions.pkl` file.

- **`Figure_3.ipynb`**  
  Reproduces Figure 3 of the main text using the baseline `all_solutions.pkl` file.

- **`Figure_4.ipynb`**  
  Reproduces Figure 4 of the main text using the baseline `all_solutions.pkl` file.

## Sensitivity tests

### `sensitivity_tests/`

Each directory contains the corresponding `venus_run_input` parameter file and, where necessary, a modified version of `run_Venus_evolution.py`.

The `all_solutions.pkl` output associated with each sensitivity test is archived on Zenodo rather than stored in this GitHub repository. After downloading the Zenodo archive, place each output file in the corresponding sensitivity-test directory to reproduce the Supplementary Figures.

| Directory | Description |
|---|---|
| `nu_1e-4/` | Core kinematic viscosity of 10⁻⁴ m² s⁻¹, used in the core–mantle friction calculation |
| `nu_10/` | Core kinematic viscosity of 10 m² s⁻¹, used in the core–mantle friction calculation |
| `Q_100/` | Constant-Q model with Q = 100 |
| `Q_30/` | Constant-Q model with Q = 30 |
| `ps_0.5bar/` | Habitable-state surface pressure of 0.5 bar |
| `ps_4bar/` | Habitable-state surface pressure of 4 bar |
| `ps_10bar/` | Habitable-state surface pressure of 10 bar |
| `SRG_280/` | Runaway-greenhouse threshold set to ASR = 280 W m⁻² |
| `SRG_320/` | Runaway-greenhouse threshold set to ASR = 320 W m⁻² |
| `qsteam_500/` | Steam-state thermal-tide amplitude q₀ = 500 Pa |
| `qsteam_4000/` | Steam-state thermal-tide amplitude q₀ = 4000 Pa |
| `OceanTides_Green_330mIT/` | Includes ocean tides in the habitable state using the tidal quality factors reported by Green et al. (2019) |
| `OceanTides_altQ/` | Includes ocean tides in the habitable state, with Q calculated from the tidal dissipation reported by Green et al. (2019) |
| `tanh_albedo_16day/` | Uses a hyperbolic-tangent albedo–rotation relationship with a critical rotation period of 16 days |
| `tanh_albedo_32day/` | Uses a hyperbolic-tangent albedo–rotation relationship with a critical rotation period of 32 days |
| `tanh_albedo_48day/` | Uses a hyperbolic-tangent albedo–rotation relationship with a critical rotation period of 48 days |

## Supplementary figures

- **`sensitivity_tests/sup_figures.py`**  
  Produces plots of habitable-state duration versus rotation period and the timing of the final habitable interval for each sensitivity test.

Usage:

```bash
python3 sup_figures.py <directory> <S_RG>
```

For example:

```bash
python3 sup_figures.py qsteam_500 300
```

This command creates the Supplementary Figure for the `qsteam_500` sensitivity test using the baseline runaway-greenhouse threshold of S<sub>RG</sub> = 300 W m⁻².

The corresponding `all_solutions.pkl` file must first be downloaded from the Zenodo archive and placed in the appropriate sensitivity-test directory.

## Reproducing the published results

There are two ways to reproduce the results:

1. **Using the archived simulation outputs:** Download the `all_solutions.pkl` files from Zenodo, place them in their corresponding directories, and run the analysis and plotting scripts in this repository.

2. **Running the simulations from scratch:** Use the model code and input files in this repository to regenerate each ensemble. Note that the published baseline ensemble contains 50,000 simulations and was run using a computing cluster.

## Data and code availability

The model source code, input files, analysis scripts, and plotting scripts are available in this GitHub repository.

The complete simulation outputs underlying the published results, including the baseline and sensitivity-test `all_solutions.pkl` files, are permanently archived on Zenodo:

**ZENODO_DOI_OR_URL**