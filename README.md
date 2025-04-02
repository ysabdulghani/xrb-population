# xrb-population

xrb-population is a Python-based project designed to study the spatial distribution of Galactic black hole low-mass X-ray binaries (BH-LMXBs) using simulated X-ray spectral modeling. The project performs extensive simulations with XSPEC to explore how various observational and intrinsic parameters affect distance estimates from soft-state (and soft-to-hard transition) spectra. In addition, it provides a framework for bias-correcting these distance estimates and comparing the corrected spatial distribution with theoretical models of the Milky Way.

## Overview

This project aims to:
- **Simulate Synthetic Spectra:** Generate fake X-ray spectra using XSPEC’s `fakeit` command with a model defined as `TBabs*(powerlaw+ezdiskbb)`.
- **Investigate Observational Biases:** Examine how parameters such as interstellar absorption (nH), photon index (Γ), disc temperature, black hole spin, mass, disc inclination, disc-to-total flux ratio, and exposure time affect the recovery of input distances.
- **Bias Correction & Analysis:** Calculate correction factors based on simulation results to mitigate systematic biases in distance estimates. The bias trends (often following an inverse-square law with distance) are then applied to correct the observed distribution.
- **Compare with Galactic Models:** Compare the bias-corrected radial distribution and galactic heights of BH-LMXBs with expected distributions (e.g., from Grimm et al. 2002), providing insights into observational selection effects and possible implications for natal kick velocities.

The simulation study supports the analysis presented in the paper:

> [Abdulghani et al. 2025, “A new independent look at the galactic black hole low-mass X-ray binary distribution”](https://arxiv.org/abs/2503.23812)

## Repository Structure

- **requirements.txt**  
  Lists the required Python packages:  
  `h5py==3.6.0`, `matplotlib==3.5.1`, `numpy==1.21.5`, `pandas==1.4.2`, and `tqdm==4.64.0`

- **xspec_simulations.py**  
  Contains the `simulation` class that:
  - Initializes simulation parameters (including instrument-specific settings for MAXI or Swift/XRT).
  - Uses XSPEC’s fakeit command to generate synthetic spectra.
  - Bins the spectrum and performs a model fit to extract flux and spectral parameters.

- **data_read.py**  
  Provides utility functions to:
  - Read date files.
  - Locate product folders based on given dates and instrument.
  - Search for saved XSPEC model files.

- **observational_effects.py**  
  Implements the simulation study to quantify observational biases by:
  - Running simulations in parallel over a grid of distances and interstellar absorption values.
  - Applying GR correction factors and normalization methods to the simulated spectra.
  - Aggregating simulation results and generating CSV tables and plots for further analysis.

- **Data Files and Database:**
  - `all_data_flat_maxi.csv` and `all_data_flat_xrt.csv`: CSV files containing simulation or observational data.
  - `results.db`: A database file used to store all simulation outcomes.

## Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/ysabdulghani/xrb-population.git
   cd xrb-population
   ```

2. **Install Dependencies:**

   It is recommended to use a virtual environment. For example:

   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows use: venv\Scripts\activate
   pip install -r requirements.txt
   ```

3. **Setup XSPEC and heasoftpy:**

   Ensure that you have XSPEC (part of the HEASoft package) installed and properly configured, as the simulations rely on PyXspec and heasoftpy.

## Usage

### Running Simulations

The main simulation study is implemented in `observational_effects.py`. You can run the script via the command line with the following arguments:

```bash
python observational_effects.py <gamma> <temp> <a> <mass> <inc> <ratio_disk_to_tot> <exposure> <instrument>
```

Where:
- `<gamma>` is the power-law photon index.
- `<temp>` is the maximum disc temperature (in keV).
- `<a>` is the black hole spin parameter (e.g., 0 or 0.998).
- `<mass>` is the black hole mass (in solar masses).
- `<inc>` is the disc inclination (in degrees).
- `<ratio_disk_to_tot>` is the disc-to-total flux ratio.
- `<exposure>` is the exposure time (in seconds).
- `<instrument>` specifies the instrument (`maxi` or `xrt`).

The script:
- Creates a temporary directory for simulation files.
- Iterates over a grid of distances and interstellar absorption (nH) values.
- Runs multiple iterations (e.g., 300 per combination) in parallel using Python’s multiprocessing.
- Saves full and reduced result tables as CSV files in the `results/<instrument>_results/` directory.

### Analyzing Results

Post-simulation, you will find CSV files summarizing:
- The fitted spectral parameters and flux estimates.
- Bias-corrected distance estimates.
- Statistical measures (e.g., red chi-squared values, fractional uncertainties).

You can then use these tables for further statistical analysis or plotting to compare with theoretical Galactic distributions.

## Scientific Context

The project is motivated by the need to understand observational biases in distance estimates to BH-LMXBs and their spatial distribution in the Milky Way. By simulating spectra over a wide parameter grid, the project investigates how systematic errors (for example, due to interstellar absorption or instrumental response) can affect the inferred distances. The bias correction approach leverages an empirical probability density function—closely resembling an inverse-square law—to de-bias the observed distribution, thereby yielding a corrected view of the Galactic population.

For more details on the methodology and scientific results, please refer to the [paper](https://arxiv.org/abs/2503.23812).

## Acknowledgements

This project makes use of:
- The Tempest High-Performance Computing Cluster at Montana State University to run the extensive multi-parameter simulations
- [XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/)
- heasoftpy for interfacing with XSPEC.
- Python libraries: NumPy, Pandas, Matplotlib, and tqdm.
- The Astropy community for numerous astronomy utilities.

For further details and discussion of the scientific results, please see the paper by [Abdulghani et al 2025](https://arxiv.org/abs/2503.23812).

