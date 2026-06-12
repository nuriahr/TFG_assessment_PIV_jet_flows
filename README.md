# TFG_assessment_PIV_jet_flows

MATLAB implementation of Event-Based Imaging Velocimetry algorithms developed for the Final Degree Project (TFG) at Universidad Carlos III de Madrid (UC3M). The algorithms are evaluated against the PaIRS ground truth on a submerged water-in-water jet experiment acquired with a Prophesee EVK4 neuromorphic camera at 200 Hz.

## Requirements
- MATLAB R2022b or later
- Image Processing Toolbox (`imdilate`, `strel`)
- Prophesee EVK4 `.hdf5` event files
- PaIRS ground truth `.mat` files

## Usage

**1. Configure paths and parameters** in `Global_Benchmark_main.m`:
```matlab
base_pairs_dir = "path/to/PaIRS_results/";
base_tifs_dir  = "path/to/TIF_images/";
datasets = {
    'Medium', "path/to/Medium_dataset.hdf5";
    'High',   "path/to/High_dataset.hdf5";
};
```

**2. Run the global benchmark:**
```matlab
run('benchmark/Global_Benchmark_main.m')
```
Outputs: RMS error table (CSV), spatial error maps (PNG),
error PDFs (PNG), and flow visualisation videos (MP4).

**3. Spectral analysis** (requires benchmark snapshots):
```matlab
run('analysis/data_point_selection.m')   % extract temporal series
run('analysis/spectrum.m')               % PSD of PaIRS signals
run('analysis/spectrum_comparison.m')    % PSD: algorithm vs. PaIRS
```

## Algorithms

| Algorithm | Description |
|---|---|
| DCC | Dense Cross-Correlation on accumulated event frames |
| FFT | FFT-based cross-correlation on Gaussian pseudo-images |
| SDC (full search) | Sparse Direct Correlation, unbounded search |
| SDC Predictor | SDC with temporal predictor warm-start |
| SDC Predictor with FFT | SDC Predictor, FFT-initialised first frame |
| SDC Predictor with Binning | SDC Predictor, binning-initialised first frame |

## Reference

Nuria Hermida Rodríguez, *Assessment of Particle Image Velocimetry using neuromorphic cameras for jet flow measurements*, TFG, Universidad Carlos III de Madrid, 2026.
