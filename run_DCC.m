function [U_val, V_val, X, Y, frame_time] = run_DCC(inFile, snap_A, snap_B, H, W, static_mask)
% RUN_DCC  Direct Cross-Correlation from an .hdf5 event file.
%
% INPUTS:
%   inFile           Path to the .hdf5 event file (string).
%   snap_A, snap_B   Frame indices (triggers) defining the two snapshots.
%   H, W             Sensor resolution [px] (e.g. 720, 1280).
%   static_mask      (optional) Logical mask (H x W). true = static-reflection
%                    pixel whose events are removed before correlation. If
%                    omitted or empty, no event filtering is applied.
%
% OUTPUTS:
%   U_val, V_val   Validated velocity field [px/frame].
%   X, Y           Grid coordinates of the interrogation window centers [px].
%   frame_time     Processing time (correlation + validation) [s].

    % Default arguments
    if nargin < 6
        static_mask = [];
    end

    % Trigger / accumulation defaults
    targetFreq = 200;
    Naccum     = 1;

    path_events = "/CD/events";
    path_trig   = "/EXT_TRIGGER/events";

    % --- 1. EVENT PARSING FROM HDF5 ---
    h5file = string(inFile);
    info       = h5info(h5file, path_events);
    NeventsTot = info.Dataspace.Size(1);

    % Final timestamp
    try
        lastEv = h5read(h5file, path_events, NeventsTot, 1);
    catch
        lastEv = h5read(h5file, path_events, [NeventsTot 1], [1 1]);
    end
    tlast = double(lastEv.t(end));

    % Hardware triggers (fall back to fixed frequency if absent)
    try
        trigs  = h5read(h5file, path_trig);
        trigOn = double(trigs.t(trigs.p == 0));
    catch
        trigOn = [];
    end

    if isempty(trigOn) || numel(trigOn) < 2
        dt     = 1e6 / targetFreq;
        trigOn = (0:dt:tlast)';
        Naccum = 1;
    else
        [~, indLastTrig] = min(abs(double(trigOn) - tlast));
        trigOn(indLastTrig+1:end) = [];
    end

    if snap_B + Naccum > numel(trigOn)
        error('Requested frames exceed the available triggers/time in the file.');
    end

    t_start_A = double(trigOn(snap_A));
    t_end_A   = double(trigOn(snap_A + Naccum));
    t_start_B = double(trigOn(snap_B));
    t_end_B   = double(trigOn(snap_B + Naccum));

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

        mask = (ev.p == 1);   % positive-polarity events only
        t_b = ev.t(mask);
        x_b = ev.x(mask);
        y_b = ev.y(mask);

        if isempty(t_b), continue; end
        if t_b(1) > t_end_B, break; end   % Frame B fully covered

        % Classify events into Frame A or Frame B by timestamp
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


    % --- 3. DIRECT CROSS-CORRELATION (DCC) ---
    tic;

    % Build binary event images
    image1 = false(H, W);
    image2 = false(H, W);

    valid_A = (x_A >= 1 & x_A <= W & y_A >= 1 & y_A <= H);
    image1(sub2ind([H, W], y_A(valid_A), x_A(valid_A))) = true;

    valid_B = (x_B >= 1 & x_B <= W & y_B >= 1 & y_B <= H);
    image2(sub2ind([H, W], y_B(valid_B), x_B(valid_B))) = true;

    % Fixed interrogation window size for all densities
    WinSize        = 64;
    OverlapPercent = 0.5;
    StepSize       = round(WinSize * (1 - OverlapPercent));

    % Adaptive minimum particle count per window
    total_events             = nnz(image1);
    avg_particles_per_window = total_events * (WinSize^2) / (H * W);
    min_particles_threshold  = max(4, round(avg_particles_per_window * 0.15));

    num_win_x = floor((W - WinSize) / StepSize) + 1;
    num_win_y = floor((H - WinSize) / StepSize) + 1;

    U_raw = nan(num_win_y, num_win_x);
    V_raw = nan(num_win_y, num_win_x);
    X     = zeros(num_win_y, num_win_x);
    Y     = zeros(num_win_y, num_win_x);

    % Full search range
    lags_y = -(WinSize-1) : (WinSize-1);
    lags_x = -(WinSize-1) : (WinSize-1);

    for wy = 1:num_win_y
        for wx = 1:num_win_x
            % Window pixel boundaries
            y_min = (wy - 1) * StepSize + 1; y_max = y_min + WinSize - 1;
            x_min = (wx - 1) * StepSize + 1; x_max = x_min + WinSize - 1;

            % Window center (for quiver plotting)
            Y(wy, wx) = (y_min + y_max) / 2;
            X(wy, wx) = (x_min + x_max) / 2;

            window1 = image1(y_min:y_max, x_min:x_max);
            window2 = image2(y_min:y_max, x_min:x_max);

            % Skip sparse background windows
            num_p = nnz(window1);
            if num_p < min_particles_threshold, continue; end

            % Correlation map over all (u, v) shifts
            corr_map = zeros(length(lags_y), length(lags_x));
            for v_idx = 1:length(lags_y)
                v = lags_y(v_idx);
                for u_idx = 1:length(lags_x)
                    u = lags_x(u_idx);
                    suma = 0;
                    for r = 1:WinSize
                        for c = 1:WinSize
                            r2 = r + v; c2 = c + u;
                            if (r2 >= 1 && r2 <= WinSize) && (c2 >= 1 && c2 <= WinSize)
                                suma = suma + (double(window1(r, c)) * double(window2(r2, c2)));
                            end
                        end
                    end
                    corr_map(v_idx, u_idx) = suma;
                end
            end

            % Adaptive correlation-peak threshold
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

                U_raw(wy, wx) = u_int + dx;
                V_raw(wy, wx) = v_int + dy;
            end
        end
    end

    time_dcc = toc;


    % --- 4. WESTERWEEL & SCARANO VALIDATION ---
    tic;
    Thr     = 2.0;
    eps_val = 0.1;
    [U_val, V_val] = WesterweelValidation(U_raw, V_raw, Thr, eps_val);
    time_val = toc;


    % --- 5. PROCESSING TIME ---
    % Event-read time excluded: it depends heavily on system load.
    frame_time = time_dcc + time_val;
end
