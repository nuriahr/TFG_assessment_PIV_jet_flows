# TFG_assessment_PIV_jet_flows

MATLAB implementation of Event-Based Imaging Velocimetry algorithms
developed as a Final Degree Project (TFG) at Universidad Carlos III de
Madrid (UC3M). Six algorithms are benchmarked against the PaIRS ground
truth on a submerged water-in-water jet acquired with a Prophesee EVK4
neuromorphic camera at 200 Hz.

## Requirements

- MATLAB R2022b or later
- Image Processing Toolbox (`imdilate`, `strel`)
- Prophesee EVK4 `.hdf5` event files
- PaIRS ground truth `.mat` files

## Scripts

### `Global_Benchmark_main.m`
Main evaluation script. For each combination of particle density and
algorithm, it computes the mean processing time, global and core RMS
error against PaIRS, time-averaged spatial error maps, flow
visualisation videos, and error probability density functions.

**Before running**, set the following variables:
```matlab
base_pairs_dir  % path to the PaIRS per-frame output folder
base_tifs_dir   % path to the folder containing the TIF image sequences
dir_error_figs  % output path for spatial error maps
dir_error_pdf   % output path for error PDF figures
dir_videos      % output path for flow videos
dir_csv         % output path for the results CSV
dir_snapshots   % output path for per-frame snapshots (spectral analysis)
datasets        % cell array with density names and .hdf5 file paths
algorithms      % cell array with algorithm names to benchmark
start_frame     % first frame to process
end_frame       % last frame to process
```

### `data_point_selection.m`
Extracts temporal velocity series at two fixed measurement points
(Shear Layer and Jet Core) from either PaIRS or algorithm snapshot
files, and saves them as `.mat` files for spectral analysis. Must be
run twice: once with `source = 'PaIRS'` and once with `source =
'<algorithm name>'`.

**Before running**, set:
```matlab
source      % 'PaIRS' or the algorithm name (e.g. 'SDC Predictor')
base_path   % path to the PaIRS results folder or snapshots folder
```

### `spectrum.m`
Computes and plots the Welch PSD of the PaIRS velocity fluctuations at
the Shear Layer and Jet Core for each density dataset. Generates three
figures per density: comparative (both points overlaid), Shear Layer
only, and Jet Core only. Requires `data_point_selection.m` to be run
first with `source = 'PaIRS'`.

**Before running**, set:
```matlab
output_dir          % output path for the PSD figures
densities           % density tags matching the data_point_selection output
U_bulk_per_density  % bulk velocity [m/s] for each density
```

### `spectrum_comparison.m`
Overlays the Welch PSD of PaIRS and a selected algorithm at the Shear
Layer and Jet Core, generating one figure per point per density.
Requires `data_point_selection.m` to be run twice beforehand (once for
PaIRS, once for the algorithm).

**Before running**, set:
```matlab
output_dir          % output path for the comparison figures
algo_name           % algorithm name matching the data_point_selection output
densities           % density tags
U_bulk_per_density  % bulk velocity [m/s] for each density
```

### Utility functions

| File | Description |
|---|---|
| `compute_bulk_velocity.m` | Computes the bulk velocity, Reynolds number, and jet centroid via axisymmetric flow-rate integration at the nozzle exit. Called automatically by `Global_Benchmark_main.m`. |
| `compute_piv_error.m` | Interpolates the evaluated velocity field onto the PaIRS reference grid and returns the global RMS error, core RMS, and spatial error map. Called automatically by `Global_Benchmark_main.m`. |
| `compute_static_mask.m` | Detects static reflection pixels in a TIF sequence by activity thresholding and morphological dilation. Called automatically by `Global_Benchmark_main.m`. |
| `WesterweelValidation.m` | Outlier detection and replacement based on the normalised median residual test (Westerweel & Scarano, 2005). Called internally by the algorithm scripts. |

### Algorithm scripts

`run_DCC.m`, `run_FFT.m`, `run_SDC.m`, `run_SDC_Predictor.m`,
`run_SDC_Predictor_FFT.m`, `run_SDC_Predictor_Binning.m`

Each script processes a pair of consecutive event frames and returns
the estimated velocity field. They are called automatically by
`Global_Benchmark_main.m` through the unified dispatcher
`run_algorithm`. They should not be run directly.

## Typical workflow
1. Run Global_Benchmark_main.m        → RMS table, error maps, PDFs, videos
2. Run data_point_selection.m         → temporal series for PaIRS (set source = 'PaIRS', enable do_save_snapshots in step 1 first)
3. Run data_point_selection.m again   → temporal series for the algorithm (set source = '<algorithm name>')
4. Run spectrum.m                     → PaIRS PSD figures
5. Run spectrum_comparison.m          → algorithm vs. PaIRS PSD figures

## Reference

Nuria Hermida Rodríguez, *Assessment of Particle Image Velocimetry using neuromorphic cameras for jet flow measurements*, TFG, Universidad Carlos III de Madrid, 2026.
