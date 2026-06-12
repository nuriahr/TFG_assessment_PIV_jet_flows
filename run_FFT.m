function [U_val, V_val, X, Y, frame_time] = run_FFT(inFile, snap_A, snap_B, H, W, static_mask)
% RUN_FFT  FFT-based cross-correlation from an .hdf5 event file.
%
%   Rasterizes events into grayscale TIF snapshots, then computes the
%   velocity field via FFT cross-correlation. Static-reflection events are
%   discarded during rasterization, so they never contribute to the images.
%
% INPUTS:
%   inFile       Path to the .hdf5 event file (string).
%   snap_A, snap_B   Frame indices (triggers) defining the two snapshots.
%   H, W         Sensor resolution [px] (e.g. 720, 1280).
%   static_mask  (optional) Logical mask (H x W). true = static-reflection
%                pixel whose events are discarded. If omitted/empty, no
%                filtering is applied.
%
% OUTPUTS:
%   U_val, V_val   Validated velocity field [px/frame].
%   X, Y           Grid coordinates of the interrogation window centers [px].
%   frame_time     Processing time (rasterization + FFT + validation) [s].

    % Default arguments
    if nargin < 6
        static_mask = [];
    end

    % Configuration
    targetFreq   = 200;
    IntPx        = 225;     % intensity added per ON event
    onlyOnEvents = true;
    useGS        = true;    % Gaussian smoothing before saving
    sigma        = 1;
    fsize        = 5;
    Naccum       = 1;
    batchSize    = 5e6;
    maxEvents    = Inf;

    snapIdxList = snap_A:snap_B;
    nSnaps      = numel(snapIdxList);

    outDir   = tempdir;          % system temp directory
    namePref = "temp_snap_";

    path_events = "/CD/events";
    path_trig   = "/EXT_TRIGGER/events";

    
    % --- 1. TIF SNAPSHOT CREATION (timed: rasterization only, not disk I/O) ---
    time_creation = 0;

    h5file = string(inFile);
    info   = h5info(h5file, path_events);
    dims   = info.Dataspace.Size;
    if numel(dims) > 1, NeventsTot = dims(1); else, NeventsTot = dims; end
    NeventsUse = min(NeventsTot, maxEvents);

    try   lastEv = h5read(h5file, path_events, NeventsUse, 1);
    catch, lastEv = h5read(h5file, path_events, [NeventsUse 1], [1 1]); end
    tlast = double(lastEv.t(end));

    try   trigs = h5read(h5file, path_trig); trigOn = trigs.t(trigs.p == 0);
    catch, trigOn = []; end

    if isempty(trigOn) || numel(trigOn) < 2
        dt     = 1e6 / targetFreq;
        trigOn = (0:dt:tlast)';
        Naccum = 1;
    else
        [~, indLastTrig] = min(abs(double(trigOn) - tlast));
        trigOn(indLastTrig+1:end) = [];
    end

    tStart_all = trigOn(snapIdxList);
    tEnd_all   = trigOn(snapIdxList + Naccum);
    clear trigOn;

    I        = zeros(H, W, "double");
    iSnap    = 1;
    startIdx = 1;

    % Extraction and rasterization loop
    while (startIdx <= NeventsUse) && (iSnap <= nSnaps)
        % Disk read (not timed)
        count = min(batchSize, NeventsUse - startIdx + 1);
        try   ev = h5read(h5file, path_events, startIdx, count);
        catch, ev = h5read(h5file, path_events, [startIdx 1], [count 1]); end
        startIdx = startIdx + count;

        % Batch processing (timed)
        tic_batch = tic;

        if onlyOnEvents
            mask = (ev.p == 1);
            t_b = ev.t(mask); x_b = ev.x(mask); y_b = ev.y(mask);
        else
            t_b = ev.t; x_b = ev.x; y_b = ev.y;
        end

        if isempty(t_b), time_creation = time_creation + toc(tic_batch); continue; end

        for k = 1:numel(t_b)
            if iSnap > nSnaps, break; end
            t = t_b(k);

            % Flush completed snapshots
            while (iSnap <= nSnaps) && (t >= tEnd_all(iSnap))
                snapNumber = snapIdxList(iSnap);
                Iout = I;
                if useGS, Iout = imgaussfilt(Iout, sigma, "FilterSize", fsize); end
                Iout = uint8(min(Iout, 255));
                fname = fullfile(outDir, sprintf("%s%05d.tif", namePref, snapNumber));
                imwrite(Iout, fname, "tif", "Compression", "lzw");

                iSnap = iSnap + 1;
                if iSnap > nSnaps, break; end
                I = zeros(H, W, "double");
            end
            if iSnap > nSnaps, break; end

            % Rasterize event into the current snapshot
            if (t > tStart_all(iSnap)) && (t < tEnd_all(iSnap))
                r = y_b(k) + 1; c = x_b(k) + 1;
                if r >= 1 && r <= H && c >= 1 && c <= W
                    % Skip events on static-reflection pixels
                    if isempty(static_mask) || ~static_mask(r, c)
                        I(r, c) = I(r, c) + IntPx;
                    end
                end
            end
        end

        time_creation = time_creation + toc(tic_batch);
    end

    % Flush remaining snapshots (timed)
    tic_cleanup = tic;
    while iSnap <= nSnaps
        snapNumber = snapIdxList(iSnap);
        Iout = I;
        if useGS, Iout = imgaussfilt(Iout, sigma, "FilterSize", fsize); end
        Iout = uint8(min(Iout, 255));
        fname = fullfile(outDir, sprintf("%s%05d.tif", namePref, snapNumber));
        imwrite(Iout, fname, "tif", "Compression", "lzw");
        iSnap = iSnap + 1;
        if iSnap <= nSnaps, I = zeros(H, W, "double"); end
    end
    time_creation = time_creation + toc(tic_cleanup);


    % --- 2. READ GENERATED TIF FILES (not timed: system-dependent) ---
    fname1 = fullfile(outDir, sprintf('%s%05d.tif', namePref, snap_A));
    fname2 = fullfile(outDir, sprintf('%s%05d.tif', namePref, snap_B));
    im1 = double(imread(fname1));
    im2 = double(imread(fname2));


    % --- 3. FFT COMPUTATION ---
    tic;

    WinSize = 64;   % fixed window size for all densities

    total_events             = nnz(im1);
    avg_particles_per_window = total_events * (WinSize^2) / (H * W);
    min_particles_threshold  = max(3, round(avg_particles_per_window * 0.15));

    Overlap = 0.5;
    step    = round(WinSize * (1 - Overlap));

    num_win_x = floor((W - WinSize) / step) + 1;
    num_win_y = floor((H - WinSize) / step) + 1;

    U_fft = nan(num_win_y, num_win_x);
    V_fft = nan(num_win_y, num_win_x);
    X     = zeros(num_win_y, num_win_x);
    Y     = zeros(num_win_y, num_win_x);

    for wy = 1:num_win_y
        for wx = 1:num_win_x
            y_min = (wy - 1) * step + 1; y_max = y_min + WinSize - 1;
            x_min = (wx - 1) * step + 1; x_max = x_min + WinSize - 1;

            Y(wy, wx) = (y_min + y_max) / 2;
            X(wy, wx) = (x_min + x_max) / 2;

            w1 = double(im1(y_min:y_max, x_min:x_max));
            w2 = double(im2(y_min:y_max, x_min:x_max));

            % Skip sparse background windows
            num_p = nnz(w1);
            if num_p < min_particles_threshold
                continue;   % U_fft/V_fft stay 0
            end

            % Zero-mean subtraction
            w1 = w1 - mean(w1(:));
            w2 = w2 - mean(w2(:));

            % Forward FFTs
            F1 = fft2(w1);
            F2 = fft2(w2);

            % Cross-power spectrum + inverse FFT
            % conj(F1).*F2 gives the correct flow direction (L<-R)
            R = real(ifft2(conj(F1) .* F2));
            R_centered = fftshift(R);

            % Primary peak
            [max_val, max_idx]   = max(R_centered(:));
            [row_peak, col_peak] = ind2sub(size(R_centered), max_idx);

            center = WinSize / 2 + 1;
            u_int  = col_peak - center;
            v_int  = row_peak - center;

            % Gaussian sub-pixel interpolation
            dx = 0; dy = 0;
            if row_peak > 1 && row_peak < WinSize && col_peak > 1 && col_peak < WinSize
                C0   = R_centered(row_peak, col_peak);
                Cm_x = R_centered(row_peak, col_peak - 1); Cp_x = R_centered(row_peak, col_peak + 1);
                Cm_y = R_centered(row_peak - 1, col_peak); Cp_y = R_centered(row_peak + 1, col_peak);

                if Cm_x > 0 && Cp_x > 0 && C0 > 0
                    dx = (log(Cm_x) - log(Cp_x)) / (2 * (log(Cm_x) - 2*log(C0) + log(Cp_x)));
                end
                if Cm_y > 0 && Cp_y > 0 && C0 > 0
                    dy = (log(Cm_y) - log(Cp_y)) / (2 * (log(Cm_y) - 2*log(C0) + log(Cp_y)));
                end
            end

            u = u_int + dx;
            v = v_int + dy;

            % Signal-to-noise ratio (mask primary peak, find secondary)
            R_mask = R_centered;
            R_mask(max(1, row_peak-2):min(WinSize, row_peak+2), ...
                   max(1, col_peak-2):min(WinSize, col_peak+2)) = 0;
            peak2_val = max(R_mask(:));
            SNR = max_val / (peak2_val + 1e-5);

            % SNR validation
            if SNR > 1.2
                U_fft(wy, wx) = u;
                V_fft(wy, wx) = v;
            end
        end
    end

    time_fft = toc;


    % --- 4. WESTERWEEL VALIDATION ---
    tic;
    Thr     = 2.0;
    eps_val = 0.1;
    [U_val, V_val] = WesterweelValidation(U_fft, V_fft, Thr, eps_val);
    time_val = toc;


    % --- 5. PROCESSING TIME + CLEANUP ---
    % Disk read time excluded (system-load dependent).
    frame_time = time_creation + time_fft + time_val;

    if isfile(fname1), delete(fname1); end
    if isfile(fname2), delete(fname2); end
end
