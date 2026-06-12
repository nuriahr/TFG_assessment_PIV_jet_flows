
% GLOBAL BENCHMARK — NEUROMORPHIC PIV
%
% Evaluates multiple particle densities and PIV algorithms against the
% PaIRS ground truth. For each (density, algorithm) pair it computes:
%   - RMS error vs. PaIRS (per frame and averaged)
%   - A spatial error map
%   - A flow visualization video
%   - The error probability density function (PDF) at three jet locations

clear; close all; clc;

set(groot, 'defaultAxesFontName',    'Times New Roman');
set(groot, 'defaultTextFontName',    'Times New Roman');
set(groot, 'defaultLegendFontName',  'Times New Roman');
set(groot, 'defaultColorbarFontName','Times New Roman');

%% --- 1. CONFIGURATION ---

% --- Directories (set before running) ---
base_pairs_dir = "path/to/PaIRS_results/";
dir_error_figs = "path/to/output/error_maps";
dir_error_pdf  = "path/to/output/error_PDFs";
dir_videos     = "path/to/output/videos";
dir_csv        = "path/to/output/csv";
dir_snapshots  = "path/to/output/snapshots";
base_tifs_dir  = "path/to/TIF_images";

if ~exist(dir_error_figs, 'dir'), mkdir(dir_error_figs); end
if ~exist(dir_videos,     'dir'), mkdir(dir_videos);     end

% --- Execution flags ---
do_save_figures = true;   % Save point-to-point error maps (.png)
do_save_videos  = true;   % Generate and save flow videos (.mp4)
do_plot_pdf     = true;   % Generate the comparative error PDF figure
do_save_csv     = true;   % Generate CSV with the computational time and the RMS
do_crop_nozzle  = false;   % true = cut the nozzle (writes NaN in U/V);
                           % false =show the complete image 720x1280
do_save_snapshots = false;  % Save per-frame U/V snapshots for spectral analysis

video_frame_step = 10;   % write one video frame every N processed frames.

% --- Datasets (set before running) ---
datasets = {
    'Sparse', "path/to/Sparse_dataset.hdf5";
    'Medium', "path/to/Medium_dataset.hdf5";
    'High',   "path/to/High_dataset.hdf5";
};

% --- Algorithms to benchmark --- 
algorithms = {'DCC', 'FFT', 'SDC (full search)', 'SDC Predictor', 'SDC Predictor with FFT', 'SDC Predictor with Binning'};

% --- Sequence and sensor configuration ---
start_frame = 2;
end_frame   = 1000;
H           = 720;     % Sensor height [px]
W           = 1280;    % Sensor width  [px]

TargetSearchRadius = 3; % [px] — used only for NaN-window reporting

% --- Physical and geometric parameters ---
cfg.x_norm_target = 1.0;              % Streamwise location to measure [x/D]
cfg.y_targets     = [0.0, 0.5, -1.0]; % Core, Shear Layer, Free Stream [y/D]
cfg.D_ext         = 160;              % Outer nozzle diameter [px]
cfg.D_int         = 120;              % Inner nozzle diameter [px]
cfg.D_ext_mm      = 30;               % Outer nozzle diameter [mm]
cfg.f_acq         = 200;              % Acquisition frequency [Hz]
cfg.nozzle_offset = 55;               % Nozzle metal offset from right edge [px]
cfg.n_calib_avg   = 10;               % Frames to average for U_bulk calibration
cfg.core_region_D = 1.0;              % Semi-width of the central region [units of D]
cfg.mask_activity_thr = 0.50;         % Fraction of frames to consider static pixel
cfg.mask_binarize_thr = 80;           % Binarization threshold [0-255]
cfg.mask_dilation_px  = 2;            % Mask dilatation radius [px]

% --- Global results table ---
global_results = table();

disp('INITIATING GLOBAL BENCHMARK...');


%% --- 2. OUTER LOOP — DENSITIES ---

for d = 1:size(datasets, 1)
    density_name = datasets{d, 1};
    hdf5_file    = datasets{d, 2};
    pairs_dir = fullfile(base_pairs_dir, density_name, 'out_PaIRS');

    fprintf('Calibrating %s for PDF analysis...\n', density_name);

    % --- Calibration: U_bulk, jet centre, measurement point indices ---
    gt_files = dir(fullfile(pairs_dir, '*.mat'));
    if isempty(gt_files)
        error('No .mat files found in: %s', pairs_dir);
    end
    fprintf('  Found %d PaIRS files in %s\n', numel(gt_files), pairs_dir);

    calib = calibrate_dataset(pairs_dir, gt_files, cfg);

    fprintf('  U_bulk = %.3f px/frame = %.3f m/s | Re = %.0f\n', ...
            calib.U_bulk, calib.U_bulk_ms, calib.Re);

    error_storage = cell(length(algorithms), 3);

    fprintf('================================================================\n');
    fprintf('>> EVALUATING DATASET: %s\n', upper(density_name));
    fprintf('================================================================\n');

    % --- Static reflection mask ---
    tifs_dir = fullfile(base_tifs_dir, ...
        sprintf('ResultsTest%s_200Hz_2026-04-28', density_name));

    if ~exist(tifs_dir, 'dir')
        error('TIF folder not found for static mask: %s', tifs_dir);
    end

    static_mask = compute_static_mask(tifs_dir, ...
        cfg.mask_activity_thr, cfg.mask_binarize_thr, cfg.mask_dilation_px);

    fprintf('  Static mask computed from: %s\n', tifs_dir);

    
    %% --- 3. INNER LOOP — ALGORITHMS ---

    for a = 1:length(algorithms)
        algo_name = algorithms{a};
        fprintf('   -> Testing algorithm: %s\n', algo_name);

        % Preallocate per-frame accumulators
        n_frames     = end_frame - start_frame;   % number of frame pairs
        rms_history  = nan(1, n_frames);
        rms_core_history = nan(1, n_frames);
        time_history = nan(1, n_frames);

        % Counters for unmeasured-window fraction in the jet region
        n_nan_jet   = 0;   % NaN windows inside the jet region (accumulated)
        n_total_jet = 0;   % total windows inside the jet region (accumulated)

        % Spatial error map accumulators
        E_map_sum   = [];
        E_map_count = [];
        X_ref_plot  = [];
        Y_ref_plot  = [];

        % Temporal predictor state (carried frame to frame, reset per algorithm)
        pred = struct('U', [], 'V', [], 'WinSize', []);

        if do_save_videos
            video_name = fullfile(dir_videos, ...
                sprintf('FlowVideo_%s_%s.mp4', density_name, algo_name));
            vid = VideoWriter(video_name, 'MPEG-4');
            vid.FrameRate = 15;
            open(vid);
            fig_vid = figure('Visible', 'off', 'Name', 'Video Generator', ...
                             'Position', [100, 100, 1000, 600]);
        end

        % --- Frame loop ---
        for snap_i = start_frame:(end_frame - 1)
            k      = snap_i - start_frame + 1;   % linear accumulator index
            snap_A = snap_i;
            snap_B = snap_i + 1;

            [U_val, V_val, X, Y, frame_time, pred] = run_algorithm( ...
                algo_name, hdf5_file, snap_A, snap_B, H, W, pred, static_mask);

            time_history(k) = frame_time;

            % Count unmeasured (NaN) windows inside the jet region
            X_alg_norm = (W - X - cfg.nozzle_offset) / cfg.D_ext;
            Y_alg_norm = (Y - calib.y_center)        / cfg.D_ext;
            jet_region = abs(Y_alg_norm) <= 1.0 & X_alg_norm >= 0 & X_alg_norm <= 6;

            n_nan_jet   = n_nan_jet   + sum(isnan(U_val(jet_region)));
            n_total_jet = n_total_jet + sum(jet_region(:));

            if do_save_snapshots
                % Build per-algorithm output folder
                snap_dir = fullfile(dir_snapshots, density_name, algo_name);
                if ~exist(snap_dir, 'dir'), mkdir(snap_dir); end

                x = X;   y = Y;   U = U_val;   V = V_val;
                snap_name = fullfile(snap_dir, sprintf('snapshot_%05d.mat', snap_i));
                save(snap_name, 'x', 'y', 'U', 'V');
            end

            % Apply global nozzle mask
            if do_crop_nozzle
                nozzle_mask = X > (max(X(:)) - cfg.nozzle_offset);
                U_val(nozzle_mask) = NaN;
                V_val(nozzle_mask) = NaN;
            end

            set(groot, 'defaultAxesFontSize',    18);
            set(groot, 'defaultTextFontSize',    18);
            set(groot, 'defaultLegendFontSize',  16);
            set(groot, 'defaultColorbarFontSize',16);

            if do_save_videos && mod(k - 1, video_frame_step) == 0
                write_video_frame(vid, fig_vid, X, Y, U_val, V_val, ...
                                  density_name, algo_name, snap_i, calib, cfg);
            end

            % --- Compare against PaIRS ground truth ---
            pairs_file = fullfile(pairs_dir, sprintf('out_%04d.mat', snap_A));

            if ~isfile(pairs_file)
                fprintf('  [WARNING] File not found: %s\n', pairs_file);
                continue;
            end
            gt = load(pairs_file);
            warning('off', 'MATLAB:interp2:NaNstrip');

            [rms_current, rms_core_current, E_map] = compute_piv_error( ...
                gt.x, gt.y, gt.U, gt.V, X, Y, U_val, V_val, calib.core_mask);
            
            rms_history(k)      = rms_current;
            rms_core_history(k) = rms_core_current;

            % Accumulate spatial error map
            if isempty(E_map_sum)
                E_map_sum   = zeros(size(E_map));
                E_map_count = zeros(size(E_map));
                X_ref_plot  = gt.x;
                Y_ref_plot  = gt.y;
            end
            valid_mask = ~isnan(E_map);
            E_map_sum(valid_mask)   = E_map_sum(valid_mask) + E_map(valid_mask);
            E_map_count(valid_mask) = E_map_count(valid_mask) + 1;

            % Capture point-wise error for the PDF (normalized by U_bulk)
            U_GT_norm  = (-double(gt.U))  / calib.U_bulk;
            U_alg_norm = (-double(U_val)) / calib.U_bulk;

            alg_x_vec = X(1, :);
            alg_y_vec = Y(:, 1);
            [~, alg_idx_x] = min(abs(alg_x_vec - calib.x_target_real));

            for p = 1:3
                [~, alg_idx_y] = min(abs(alg_y_vec - calib.y_target_real(p)));
                err = U_alg_norm(alg_idx_y, alg_idx_x) - ...
                      U_GT_norm(calib.p_idx_y(p), calib.p_idx_x);
                error_storage{a, p}(end+1) = err;
            end

            if mod(k, 10) == 0 || snap_i == end_frame - 1
                fprintf('      %s | %s: frame %d/%d\n', ...
                        density_name, algo_name, k, n_frames);
            end

        end % frame loop

        fprintf('  Errors captured for PDF — Core:%d, Shear:%d, FreeStream:%d\n', ...
                numel(error_storage{a,1}), numel(error_storage{a,2}), numel(error_storage{a,3}));

        if do_save_videos
            close(vid); close(fig_vid);
            fprintf('   -> Video saved to: %s\n', video_name);
        end

        % --- Finalize algorithm results ---
        if any(~isnan(rms_history))
            rms_mean = mean(rms_history, 'omitnan');
            rms_core_mean = mean(rms_core_history, 'omitnan');
            rms_mean_norm  = rms_mean      / calib.U_bulk;
            rms_core_norm  = rms_core_mean / calib.U_bulk;

            % Report unmeasured-window fraction in the jet region
            if n_total_jet > 0
                nan_frac_jet = 100 * n_nan_jet / n_total_jet;
                fprintf('   [%s] Unmeasured windows in jet region: %.1f%% (radius = %d px)\n', ...
                        algo_name, nan_frac_jet, TargetSearchRadius); 
            end

            % Skip first frame for timing (cold start)
            valid_times = time_history(~isnan(time_history));
            if numel(valid_times) > 1
                time_mean = mean(valid_times(2:end));
            else
                time_mean = valid_times(1);
            end

            global_results = [global_results; ...
                {density_name, algo_name, time_mean, ...
                 rms_mean, rms_mean_norm, rms_core_mean, rms_core_norm}];

            % Averaged spatial error map (normalized by U_bulk)
            E_map_mean = E_map_sum ./ max(E_map_count, 1);
            E_map_mean(E_map_count == 0) = NaN;
            E_map_norm = E_map_mean / calib.U_bulk;
 
            if do_save_figures
                % Normalized coordinates (x/D, y/D)
                X_norm = (max(X_ref_plot(:)) - X_ref_plot - cfg.nozzle_offset) / cfg.D_ext;
                Y_norm = (Y_ref_plot - calib.y_center) / cfg.D_ext;

                set(groot, 'defaultAxesFontSize',    20);
                set(groot, 'defaultTextFontSize',    20);
                set(groot, 'defaultLegendFontSize',  14);
                set(groot, 'defaultColorbarFontSize',14);
     
                fig = figure('Visible', 'off', 'Name', [algo_name, ' - ', density_name]);
                contourf(X_norm, Y_norm, E_map_norm, 20, 'LineColor', 'none');
                set(gca, 'YDir', 'reverse');
                colormap(jet);
                set(gca, 'CLim', [0, 0.5]);
                c = colorbar; 
                c.Label.String = 'Average Error |\Delta U| / U_{bulk}';
                c.Label.FontSize = 16;
                c.FontSize = 16;
                title(sprintf('RMS = %.3f | %.2f px', rms_mean / calib.U_bulk, rms_mean));
                xlabel('x/D');
                ylabel('y/D');
                axis equal; 
                axis tight;
                xlim([0, 7]);
                % ylim([min(Y_norm(:)), max(Y_norm(:))]);
                ylim([-1.9, 1.9]);
                set(gca, 'YTick', [-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5]);
                fig_name = fullfile(dir_error_figs, ...
                    sprintf('ErrorMap_%s_%s.png', density_name, algo_name));
                saveas(fig, fig_name); close(fig);
                fprintf('  -> Figure saved to: %s\n', fig_name);
            end
        end
    end % algorithm loop

    %% --- 4. COMPARATIVE ERROR PDF FIGURE ---
    set(groot, 'defaultAxesFontSize',    18);
    set(groot, 'defaultTextFontSize',    18);
    set(groot, 'defaultLegendFontSize',  16);
    set(groot, 'defaultColorbarFontSize',16);

    if do_plot_pdf
        plot_error_pdf(error_storage, algorithms, density_name, ...
                       cfg.x_norm_target, dir_error_pdf, do_save_figures, calib.U_bulk);
    end

end % density loop

%% --- 5. EXPORT RESULTS TABLE ---
fprintf('================================================================\n');
fprintf('                  GLOBAL BENCHMARK SUMMARY                       \n');
fprintf('================================================================\n');

if ~isempty(global_results)
    global_results.Properties.VariableNames = ...
        {'Density', 'Algorithm', 'Mean_Time_s', ...
         'RMS_px', 'RMS_norm', 'RMS_core_px', 'RMS_core_norm'};
    disp(global_results);

    if do_save_csv
        csv_path = fullfile(dir_csv, 'benchmark_results.csv');
        writetable(global_results, csv_path);
        fprintf('-> Results saved to: %s\n', csv_path);
    end
else
    disp('No results were obtained.');
end


%%  LOCAL FUNCTIONS

function [U_val, V_val, X, Y, frame_time, pred] = run_algorithm( ...
    algo_name, hdf5_file, snap_A, snap_B, H, W, pred, static_mask)
% RUN_ALGORITHM  Unified dispatcher for all PIV algorithms.
%   Predictor-based algorithms read and update the predictor state struct
%   'pred' (fields .U, .V, .WinSize). Non-predictor algorithms ignore it.

    switch algo_name
        case 'DCC'
            [U_val, V_val, X, Y, frame_time] = ...
                run_DCC(hdf5_file, snap_A, snap_B, H, W, static_mask);

        case 'FFT'
            [U_val, V_val, X, Y, frame_time] = ...
                run_FFT(hdf5_file, snap_A, snap_B, H, W, static_mask);

        case 'SDC (full search)'
            [U_val, V_val, X, Y, frame_time] = ...
                run_SDC(hdf5_file, snap_A, snap_B, H, W, static_mask);

        case 'SDC Predictor'
            [U_val, V_val, X, Y, pred.U, pred.V, frame_time] = ...
                run_SDC_Predictor(hdf5_file, snap_A, snap_B, H, W, ...
                    pred.U, pred.V, static_mask);

        case 'SDC Predictor with FFT'
            [U_val, V_val, X, Y, pred.U, pred.V, frame_time] = ...
                run_SDC_Predictor_FFT(hdf5_file, snap_A, snap_B, H, W, ...
                    pred.U, pred.V, static_mask);

        case 'SDC Predictor with Binning'
            [U_val, V_val, X, Y, pred.U, pred.V, frame_time] = ...
                run_SDC_Predictor_Binning(hdf5_file, snap_A, snap_B, H, W, ...
                    pred.U, pred.V, static_mask);

        otherwise
            error('run_algorithm: unknown algorithm "%s"', algo_name);
    end
end

function write_video_frame(vid, fig_vid, X, Y, U_val, V_val, ...
                           density_name, algo_name, snap_i, calib, cfg)
% WRITE_VIDEO_FRAME  Renders one normalized velocity-magnitude frame with
%   vectors and writes the video.

    clf(fig_vid);
    ax = axes(fig_vid);

    % Normalized coordinates (x/D, y/D), nozzle at x/D = 0
    X_norm = (max(X(:)) - X - cfg.nozzle_offset) / cfg.D_ext;
    Y_norm = (Y - calib.y_center) / cfg.D_ext;

    % Normalized velocity components
    U_norm = -U_val / calib.U_bulk;
    V_norm =  V_val / calib.U_bulk;

    % Normalized velocity magnitude
    vel_mag = sqrt(U_norm.^2 + V_norm.^2);

    pcolor(ax, X_norm, Y_norm, vel_mag);
    shading(ax, 'interp');
    hold(ax, 'on');
    quiver(ax, X_norm, Y_norm, U_norm, V_norm, 1.5, 'Color', 'k', 'LineWidth', 0.8);
    hold(ax, 'off');

    % Y axis still points downward (image convention); X grows with the flow
    set(ax, 'YDir', 'reverse');
    axis(ax, 'equal'); axis(ax, 'tight');
    xlim(ax, [0, max(X_norm(:))]);
    ylim(ax, [min(Y_norm(:)), max(Y_norm(:))]);

    colormap(ax, turbo);
    c = colorbar(ax);
    c.Label.String = '|U| / U_{bulk}';
    clim(ax, [0, 1.5]);
    title(ax, sprintf('%s: %s | Frame %d', density_name, algo_name, snap_i));
    xlabel(ax, 'x/D'); ylabel(ax, 'y/D'); grid(ax, 'on');

    writeVideo(vid, getframe(fig_vid));
end


function calib = calibrate_dataset(pairs_dir, gt_files, cfg)
% CALIBRATE_DATASET  Computes U_bulk and the measurement-point indices from
%   the averaged PaIRS field. Averaging several frames makes the centroid
%   and bulk velocity robust against single-frame noise.

    % Load the statistically converged PaIRS average field (out.mat)
    avg_file = fullfile(pairs_dir, 'out.mat');
    if isfile(avg_file)
        avg_data = load(avg_file, 'x', 'y', 'U');
        x      = double(avg_data.x);
        y      = double(avg_data.y);
        U_mean = double(avg_data.U);
        fprintf('  [calibrate_dataset] Loaded converged average: out.mat\n');
    else
        % Fallback: average first n_calib_avg instantaneous frames
        fprintf('  [calibrate_dataset] out.mat not found, averaging %d frames\n', cfg.n_calib_avg);
        load(fullfile(pairs_dir, gt_files(1).name), 'x', 'y', 'U');
        n_avg  = min(cfg.n_calib_avg, numel(gt_files));
        U_sum  = zeros(size(U));
        for ka = 1:n_avg
            tmp   = load(fullfile(pairs_dir, gt_files(ka).name), 'U');
            U_sum = U_sum + double(tmp.U);
        end
        U_mean = U_sum / n_avg;
    end

    % Coordinate transform: X=0 at nozzle exit, U positive (flow L→R)
    x_shifted = max(x(:)) - x - cfg.nozzle_offset;
    U_flipped = -U_mean;

    % Bulk velocity calculation
    [calib.U_bulk, calib.U_bulk_ms, calib.Re, bulk_diag] = compute_bulk_velocity( ...
        x_shifted, y, U_flipped, cfg.D_int, cfg.D_ext, cfg.D_ext_mm, cfg.f_acq);
    calib.y_center = bulk_diag.y_center_px;

    % Locate the three measurement points on the PaIRS grid
    x_n = x_shifted / cfg.D_ext;
    y_n = (y - calib.y_center) / cfg.D_ext;

    [~, calib.p_idx_x] = min(abs(x_n(1,:) - cfg.x_norm_target));
    calib.p_idx_y = zeros(1, 3);
    for p = 1:3
        [~, calib.p_idx_y(p)] = min(abs(y_n(:,1) - cfg.y_targets(p)));
    end

    % Real physical coordinates [px]
    calib.x_target_real = x(1, calib.p_idx_x);
    calib.y_target_real = arrayfun(@(p) y(calib.p_idx_y(p), 1), 1:3);

    % Mask of the center region (jet core) over PaIRS grid
    gt_x_n = x_shifted / cfg.D_ext;
    gt_y_n = (y - calib.y_center) / cfg.D_ext;
    calib.core_mask = abs(gt_y_n) <= cfg.core_region_D & ...
                      abs(gt_x_n - cfg.x_norm_target) <= 0.5;
end


function plot_error_pdf(error_storage, algorithms, density_name, ...
                        x_norm_target, dir_error_figs, do_save, U_bulk)
% PLOT_ERROR_PDF  Builds the 3-panel comparative error PDF figure
%   (Jet Core, Shear Layer, Free Stream) for all algorithms.

    f_pdf = figure('Name', ['PDF Error - ', density_name], ...
                   'Color', 'w', 'Position', [100 100 1300 450]);

    num_algos   = length(algorithms);
    colors      = turbo(max(num_algos, 1));
    line_styles = {'-', '--', ':', '-.'};
    markers     = {'o', 's', 'd', '^', 'v', 'p', 'h'};
    point_titles = {'Jet Core (y/D=0.0)', 'Shear Layer (y/D=0.5)', ...
                    'Free Stream (y/D=-1.0)'};

    for p = 1:3
        ax_p = subplot(1, 3, p, 'Parent', f_pdf);

        pos = get(ax_p, 'Position');
        pos(2) = pos(2) + 0.20;
        pos(4) = pos(4) - 0.34;
        set(ax_p, 'Position', pos);

        cla(ax_p);
        hold(ax_p, 'on');

        legend_entries = {};
        pdf_written    = false;

        for a = 1:num_algos
            point_data = error_storage{a, p};

            if ~isempty(point_data) && numel(point_data) > 10 && ~all(isnan(point_data))
                try
                    sigma = std(point_data, 'omitnan');
                    [f_vals, x_vals] = ksdensity(point_data);

                    plot(ax_p, x_vals, f_vals, ...
                        'Color',           colors(a, :), ...
                        'LineStyle',       line_styles{mod(a-1, numel(line_styles)) + 1}, ...
                        'LineWidth',       1.5, ...
                        'Marker',          markers{mod(a-1, numel(markers)) + 1}, ...
                        'MarkerSize',      6, ...
                        'MarkerFaceColor', 'none', ...
                        'MarkerIndices',   round(linspace(1, length(x_vals), 15)));

                    legend_entries{end+1} = sprintf('%s (std:%.3f)', algorithms{a}, sigma);
                    pdf_written = true;
                catch
                    warning('ksdensity failed for point %d, algorithm %s', p, algorithms{a});
                end
            end
        end

       if pdf_written
            xline(ax_p, 0, 'k--', 'HandleVisibility', 'off');
            grid(ax_p, 'on');

            t = title(ax_p, point_titles{p});
            t.Units = 'normalized';
            t.Position(2) = 1.20;

            xlabel(ax_p, '\Delta U / U_{bulk}');
            if p == 1, ylabel(ax_p, 'Probability Density'); end
            colormap(ax_p, 'default');
        end
        hold(ax_p, 'off');

        % Secondary x-axis in pixels
        if pdf_written
            xlims_norm = xlim(ax_p);
        
            ax_top = axes('Position',        get(ax_p, 'Position'), ...
                          'XAxisLocation',   'top',                 ...
                          'YAxisLocation',   'right',               ...
                          'Color',           'none',                ...
                          'YTick',           [],                    ...
                          'YColor',          'none',                ...
                          'FontName',        'Times New Roman',     ...
                          'FontSize',        18,                    ...
                          'Parent',          f_pdf);
        
            ax_top.XLim = xlims_norm * U_bulk;
            xlabel(ax_top, '\Delta U (px)');
        
            ax_top.XLabel.Units = 'normalized';
            ax_top.XLabel.Position(2) = 1.10;
        
            uistack(ax_top, 'top');
        end

    end   % for p = 1:3

    lgd = legend(ax_p, legend_entries, 'Interpreter', 'none', ...
             'Orientation', 'horizontal', 'NumColumns', 2, ...
             'FontSize', 16);
    lgd.Units = 'normalized';
    drawnow;
    lgd.Position(1) = 0.5 - lgd.Position(3)/2;
    lgd.Position(2) = 0.01;

    if do_save
        pdf_fig_name = fullfile(dir_error_figs, sprintf('PDF_Error_%s.png', density_name));
        saveas(f_pdf, pdf_fig_name);
        fprintf('  -> PDF figure saved to: %s\n', pdf_fig_name);
    end
end
