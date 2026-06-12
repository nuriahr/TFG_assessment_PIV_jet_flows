% DATA_POINT_SELECTION
% Extracts temporal velocity series at the Shear Layer and Jet Core and
% saves them as .mat files for downstream spectral analysis.
% Run once with source = 'PaIRS' and once per algorithm to compare.

clear; clc; close all;

%% 1. PARAMETERS
D_ext_px    = 160;      % [px]
D_ext_mm    = 30;       % [mm]
x_offset_px = 55;       % [px] nozzle offset from right sensor edge
fs          = 200;      % [Hz]
px_to_ms    = (D_ext_mm * 1e-3) / D_ext_px / (1/fs);   % [px/frame] -> [m/s]

densities = {'Medium', 'High'};

% Source: 'PaIRS' reads PaIRS per-frame outputs; algorithm name reads
% snapshots saved by Global_Benchmark_main with do_save_snapshots = true.
source = 'PaIRS';

if strcmp(source, 'PaIRS')
    base_path  = 'path/to/Results_PaIRS';
    subfolder  = 'out_PaIRS';
    source_tag = 'PaIRS';
else
    base_path  = 'path/to/Snapshots';
    subfolder  = source;
    source_tag = matlab.lang.makeValidName(source);
end

% Measurement points [x/D, y/D]
x_SL_target_norm   = 3.0;   y_SL_target_norm   =  0.5;   % Shear Layer
x_Core_target_norm = 2.0;   y_Core_target_norm =  0.0;   % Jet Core


%% 2. MAIN LOOP OVER DENSITIES
for d = 1:length(densities)
    type       = densities{d};
    type_lower = lower(type);
    fprintf('\n=== Processing density: %s ===\n', type);

    folder_data = fullfile(base_path, type, subfolder);
    files = dir(fullfile(folder_data, '*.mat'));
    if isempty(files)
        warning('No .mat files found in: %s — skipping.', folder_data);
        continue;
    end

    % --- Jet centroid from average of first 10 frames ---
    load(fullfile(folder_data, files(1).name), 'x', 'y', 'U');
    n_avg = min(10, length(files));
    U_sum = zeros(size(U));
    for ka = 1:n_avg
        tmp   = load(fullfile(folder_data, files(ka).name), 'U');
        U_sum = U_sum + tmp.U;
    end
    U_mean_flipped = -(U_sum / n_avg);   % flip sign: raw flow is right-to-left

    x_mid = max(x(1,:)) - x_offset_px - 1.0 * D_ext_px;
    [~, idx_x_mid]    = min(abs(x(1,:) - x_mid));
    [~, idx_y_center] = max(U_mean_flipped(:, idx_x_mid));
    y_center_px       = y(idx_y_center, idx_x_mid);

    % --- Grid indices of the measurement points ---
    x_target_SL   = max(x(1,:)) - x_offset_px - x_SL_target_norm   * D_ext_px;
    x_target_Core = max(x(1,:)) - x_offset_px - x_Core_target_norm * D_ext_px;

    [~, idx_x_SL]   = min(abs(x(1,:) - x_target_SL));
    [~, idx_y_SL]   = min(abs(y(:,1) - (y_center_px + y_SL_target_norm   * D_ext_px)));
    [~, idx_x_Core] = min(abs(x(1,:) - x_target_Core));
    [~, idx_y_Core] = min(abs(y(:,1) - (y_center_px + y_Core_target_norm * D_ext_px)));

    fprintf('  Centroid:    y = %.1f px\n',              y_center_px);
    fprintf('  Shear Layer: x = %.1f px, y = %.1f px\n', x(1,idx_x_SL),   y(idx_y_SL,1));
    fprintf('  Jet Core:    x = %.1f px, y = %.1f px\n', x(1,idx_x_Core), y(idx_y_Core,1));

    % --- Temporal series extraction (U negated; V keeps sign) ---
    num_frames  = length(files);
    u_time_SL   = zeros(1, num_frames);
    v_time_SL   = zeros(1, num_frames);
    u_time_Core = zeros(1, num_frames);
    v_time_Core = zeros(1, num_frames);

    for i = 1:num_frames
        data           = load(fullfile(folder_data, files(i).name), 'U', 'V');
        u_time_SL(i)   = -data.U(idx_y_SL,   idx_x_SL)   * px_to_ms;
        v_time_SL(i)   =  data.V(idx_y_SL,   idx_x_SL)   * px_to_ms;
        u_time_Core(i) = -data.U(idx_y_Core, idx_x_Core)  * px_to_ms;
        v_time_Core(i) =  data.V(idx_y_Core, idx_x_Core)  * px_to_ms;
    end

    % --- Fill NaN frames by linear interpolation ---
    n_nan = sum(isnan(u_time_SL) | isnan(v_time_SL) | ...
                isnan(u_time_Core) | isnan(v_time_Core));
    if n_nan > 0
        fprintf('  NaN frames: %d / %d — filling by linear interpolation\n', ...
                n_nan, num_frames);
        u_time_SL   = fillmissing(u_time_SL,   'linear');
        v_time_SL   = fillmissing(v_time_SL,   'linear');
        u_time_Core = fillmissing(u_time_Core, 'linear');
        v_time_Core = fillmissing(v_time_Core, 'linear');
    end

    % --- Save ---
    name_SL   = sprintf('data_spectrum_ShearLayer_%s_%s.mat', type_lower, source_tag);
    name_Core = sprintf('data_spectrum_Core_%s_%s.mat',       type_lower, source_tag);

    save(name_SL,   'u_time_SL',   'v_time_SL',   'fs', ...
                    'x_SL_target_norm',   'y_SL_target_norm');
    save(name_Core, 'u_time_Core', 'v_time_Core', 'fs', ...
                    'x_Core_target_norm', 'y_Core_target_norm');

    fprintf('  Saved: %s\n  Saved: %s\n', name_SL, name_Core);
end

disp('--- Data point selection complete ---');