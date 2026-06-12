function [U_val, V_val, X, Y, U_pred_out, V_pred_out, frame_time] = ...
    run_SDC_Predictor_Binning(inFile, snap_A, snap_B, H, W, U_pred_in, V_pred_in, static_mask)
% RUN_SDC_PREDICTOR_BINNING  SDC with a temporal predictor and a fast
%   2x2-binned warm start.
%
%   On the first frame the search is performed on a 2x2-binned image, which
%   halves the search space and speeds up the cold start; the resulting
%   displacement is scaled back by 2 to full resolution. Subsequent frames
%   use the temporal predictor with a targeted search at full resolution.
%
% INPUTS:
%   inFile                 Path to the .hdf5 event file (string).
%   snap_A, snap_B         Frame indices (triggers) defining the two snapshots.
%   H, W                   Sensor resolution [px] (e.g. 720, 1280).
%   U_pred_in, V_pred_in   Predictor field from the previous frame
%                          (empty [] triggers a binned warm start).
%   WinSize_in             Window size propagated from the previous frame.
%   static_mask            (optional) Logical mask (H x W). true = static-
%                          reflection pixel whose events are removed before
%                          correlation. If omitted/empty, no filtering applied.
%
% OUTPUTS:
%   U_val, V_val             Validated velocity field [px/frame].
%   X, Y                     Grid coordinates of the window centers [px].
%   U_pred_out, V_pred_out   Predictor field for the next frame.
%   WinSize_out              Window size to propagate to the next frame.
%   frame_time               Processing time (correlation + validation) [s].

    % Default arguments
    if nargin < 8
        static_mask = [];
    end

    % Trigger / accumulation defaults
    targetFreq         = 200;
    Naccum             = 1;

    TargetSearchRadius = 3;    % search radius around the temporal predictor
    WarmStartRadius    = 20;   % full-search radius on the warm start [px]

    path_events = "/CD/events";
    path_trig   = "/EXT_TRIGGER/events";


    % --- 1. EVENT PARSING FROM HDF5 ---
    h5file = string(inFile);
    info       = h5info(h5file, path_events);
    NeventsTot = info.Dataspace.Size(1);

    try   lastEv = h5read(h5file, path_events, NeventsTot, 1);
    catch, lastEv = h5read(h5file, path_events, [NeventsTot 1], [1 1]); end
    tlast = double(lastEv.t(end));

    try   trigs = h5read(h5file, path_trig); trigOn = double(trigs.t(trigs.p == 0));
    catch, trigOn = []; end

    if isempty(trigOn) || numel(trigOn) < 2
        dt     = 1e6 / targetFreq;
        trigOn = (0:dt:tlast)';
        Naccum = 1;
    else
        [~, indLastTrig] = min(abs(double(trigOn) - tlast));
        trigOn(indLastTrig+1:end) = [];
    end

    t_start_A = double(trigOn(snap_A)); t_end_A = double(trigOn(snap_A + Naccum));
    t_start_B = double(trigOn(snap_B)); t_end_B = double(trigOn(snap_B + Naccum));

    % Preallocate coordinate lists (sparse event representation)
    max_est = 2000000;
    x_A = zeros(max_est, 1); y_A = zeros(max_est, 1); count_A = 0;
    x_B = zeros(max_est, 1); y_B = zeros(max_est, 1); count_B = 0;

    batchSize = 5e6;
    startIdx  = 1;
    while startIdx <= NeventsTot
        count = min(batchSize, NeventsTot - startIdx + 1);
        ev = h5read(h5file, path_events, startIdx, count);
        startIdx = startIdx + count;

        mask = (ev.p == 1);
        t_b = ev.t(mask); x_b = ev.x(mask); y_b = ev.y(mask);
        if isempty(t_b), continue; end
        if double(t_b(1)) > t_end_B, break; end

        for k = 1:numel(t_b)
            t = double(t_b(k));
            if (t > t_start_A) && (t <= t_end_A)
                count_A = count_A + 1; x_A(count_A) = x_b(k); y_A(count_A) = y_b(k);
            elseif (t > t_start_B) && (t <= t_end_B)
                count_B = count_B + 1; x_B(count_B) = x_b(k); y_B(count_B) = y_b(k);
            end
        end
    end

    % Trim to actual size and remove duplicate pixels
    coords_A = unique([x_A(1:count_A), y_A(1:count_A)], 'rows');
    if isempty(coords_A), x_A = []; y_A = []; else, x_A = coords_A(:,1); y_A = coords_A(:,2); end

    coords_B = unique([x_B(1:count_B), y_B(1:count_B)], 'rows');
    if isempty(coords_B), x_B = []; y_B = []; else, x_B = coords_B(:,1); y_B = coords_B(:,2); end


    % --- 2. STATIC-REFLECTION FILTERING ---
    % Remove events that fall on static-reflection pixels. 
    % Events use 0-based coordinates; add 1 to index the MATLAB mask.
    if ~isempty(static_mask)
        if ~isempty(x_A)
            keep_A = ~static_mask(sub2ind([H, W], y_A + 1, x_A + 1));
            x_A = x_A(keep_A);  y_A = y_A(keep_A);
        end
        if ~isempty(x_B)
            keep_B = ~static_mask(sub2ind([H, W], y_B + 1, x_B + 1));
            x_B = x_B(keep_B);  y_B = y_B(keep_B);
        end
    end


    % --- 3. BINNING + SPARSE DIRECT CORRELATION + TEMPORAL PREDICTOR ---
    tic_sdc = tic;

    % Build binary event images
    image1 = false(H, W);
    image2 = false(H, W);

    valid_A = (x_A >= 1 & x_A <= W & y_A >= 1 & y_A <= H);
    image1(sub2ind([H, W], y_A(valid_A), x_A(valid_A))) = true;

    valid_B = (x_B >= 1 & x_B <= W & y_B >= 1 & y_B <= H);
    image2(sub2ind([H, W], y_B(valid_B), x_B(valid_B))) = true;

    % First frame: binned warm start
    is_warm_start = isempty(U_pred_in);

    % Boundary condition (spatial extrapolation)
    % Flow runs right-to-left: particles enter through the last column,
    % which has no history. Copy the velocity from the adjacent column.
    if ~is_warm_start && size(U_pred_in, 2) > 1
        U_pred_in(:, end) = U_pred_in(:, end-1);
        V_pred_in(:, end) = V_pred_in(:, end-1);
    end

    total_events_current = nnz(image1);
     WinSize     = 64;

    % Window size and warm-start binning setup
    if is_warm_start

        % Apply 2x2 binning to halve the warm-start search space
        WinSize_bin = max(16, round(WinSize / 2));
        WinSize_bin = WinSize_bin + mod(WinSize_bin, 2);

        im1_bin = image1(1:2:end, 1:2:end) | image1(2:2:end, 1:2:end) | ...
                  image1(1:2:end, 2:2:end) | image1(2:2:end, 2:2:end);
        im2_bin = image2(1:2:end, 1:2:end) | image2(2:2:end, 1:2:end) | ...
                  image2(1:2:end, 2:2:end) | image2(2:2:end, 2:2:end);

        % Full search range in binned space
        % lags_y = -(WinSize_bin-1) : (WinSize_bin-1);
        % lags_x = -(WinSize_bin-1) : (WinSize_bin-1);

        % Bounded full search range in binned space
        r_bin  = min(round(WarmStartRadius/2), WinSize_bin-1);
        lags_y = -r_bin : r_bin;
        lags_x = -r_bin : r_bin;

    end

    OverlapPercent = 0.5;
    StepSize       = round(WinSize * (1 - OverlapPercent));

    num_win_x = floor((W - WinSize) / StepSize) + 1;
    num_win_y = floor((H - WinSize) / StepSize) + 1;

    avg_particles_per_window = total_events_current * (WinSize^2) / (H * W);
    min_particles_threshold  = max(4, round(avg_particles_per_window * 0.15));

    U_raw = nan(num_win_y, num_win_x);
    V_raw = nan(num_win_y, num_win_x);
    X     = zeros(num_win_y, num_win_x);
    Y     = zeros(num_win_y, num_win_x);

    for wy = 1:num_win_y
        for wx = 1:num_win_x
            % Full-resolution window boundaries
            y_min = (wy - 1) * StepSize + 1; y_max = y_min + WinSize - 1;
            x_min = (wx - 1) * StepSize + 1; x_max = x_min + WinSize - 1;

            Y(wy, wx) = (y_min + y_max) / 2;
            X(wy, wx) = (x_min + x_max) / 2;

            if is_warm_start
                % Map boundaries to binned space
                y_min_bin = floor((y_min-1)/2) + 1; y_max_bin = y_min_bin + WinSize_bin - 1;
                x_min_bin = floor((x_min-1)/2) + 1; x_max_bin = x_min_bin + WinSize_bin - 1;

                % Clamp to binned-image bounds
                y_max_bin = min(y_max_bin, size(im1_bin, 1));
                x_max_bin = min(x_max_bin, size(im1_bin, 2));

                window1 = im1_bin(y_min_bin:y_max_bin, x_min_bin:x_max_bin);
                window2 = im2_bin(y_min_bin:y_max_bin, x_min_bin:x_max_bin);
                WinSize_current = WinSize_bin;
            else
                % Full-resolution predictor with targeted search
                window1 = image1(y_min:y_max, x_min:x_max);
                window2 = image2(y_min:y_max, x_min:x_max);
                WinSize_current = WinSize;
                
                if isnan(U_pred_in(wy, wx)) || isnan(V_pred_in(wy, wx))
                    % Predictor missing here: use warm start radius
                    r      = min(WarmStartRadius, WinSize-1);   % full-res (±20 px) because no warm start in this case
                    lags_y = -r : r;
                    lags_x = -r : r;

                else
                    u_pred = round(U_pred_in(wy, wx));
                    v_pred = round(V_pred_in(wy, wx));
                    lags_x = (u_pred - TargetSearchRadius) : (u_pred + TargetSearchRadius);
                    lags_y = (v_pred - TargetSearchRadius) : (v_pred + TargetSearchRadius);
                end
            end

            [p_r, p_c] = find(window1);
            num_p = numel(p_r);
            if num_p < min_particles_threshold, continue; end

            % SDC core: iterate over events only, not pixels
            corr_map = zeros(length(lags_y), length(lags_x));
            for v_idx = 1:length(lags_y)
                v = lags_y(v_idx);
                for u_idx = 1:length(lags_x)
                    u = lags_x(u_idx);
                    suma = 0;
                    for k = 1:num_p
                        r2 = p_r(k) + v; c2 = p_c(k) + u;
                        if (r2 >= 1 && r2 <= WinSize_current) && (c2 >= 1 && c2 <= WinSize_current)
                            if window2(r2, c2), suma = suma + 1; end
                        end
                    end
                    corr_map(v_idx, u_idx) = suma;
                end
            end

            min_corr_peak = max(3, round(num_p * 0.15));
            [max_val, max_idx] = max(corr_map(:));

            if max_val >= min_corr_peak
                [row_peak, col_peak] = ind2sub(size(corr_map), max_idx);

                u_int = lags_x(col_peak);
                v_int = lags_y(row_peak);

                % Gaussian sub-pixel interpolation
                dx = 0; dy = 0;
                if row_peak > 1 && row_peak < size(corr_map, 1) && ...
                   col_peak > 1 && col_peak < size(corr_map, 2)
                    C0   = corr_map(row_peak, col_peak);
                    Cm_x = corr_map(row_peak, col_peak - 1); Cp_x = corr_map(row_peak, col_peak + 1);
                    Cm_y = corr_map(row_peak - 1, col_peak); Cp_y = corr_map(row_peak + 1, col_peak);
                    if Cm_x > 0 && Cp_x > 0 && C0 > 0
                        dx = (log(Cm_x) - log(Cp_x)) / (2 * (log(Cm_x) - 2*log(C0) + log(Cp_x)));
                    end
                    if Cm_y > 0 && Cp_y > 0 && C0 > 0
                        dy = (log(Cm_y) - log(Cp_y)) / (2 * (log(Cm_y) - 2*log(C0) + log(Cp_y)));
                    end
                end

                % On the warm start the result is in binned space:
                % scale displacement back to full resolution (x2).
                if is_warm_start
                    U_raw(wy, wx) = (u_int + dx) * 2;
                    V_raw(wy, wx) = (v_int + dy) * 2;
                else
                    U_raw(wy, wx) = u_int + dx;
                    V_raw(wy, wx) = v_int + dy;
                end
            end
        end
    end

    time_sdc = toc(tic_sdc);


    % --- 4. WESTERWEEL VALIDATION + PREDICTOR UPDATE ---
    tic_val = tic;
    Thr     = 2.0;
    eps_val = 0.1;
    [U_val, V_val] = WesterweelValidation(U_raw, V_raw, Thr, eps_val);

    % --- Nozzle mask ---
    % Force NaN in the nozzle region so the static metal does not
    % contaminate the predictor passed to the next frame.
    offset_boquera = 55;
    nozzle_mask = X > (W - offset_boquera);
    U_val(nozzle_mask) = NaN;
    V_val(nozzle_mask) = NaN;

    % Propagate the validated field as the next frame's predictor
    U_pred_out = U_val;
    V_pred_out = V_val;

    time_val = toc(tic_val);


    % --- 5. PROCESSING TIME ---
    % Event-read time excluded: it depends heavily on system load.
    frame_time = time_sdc + time_val;
end
