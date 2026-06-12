function [U_bulk, U_bulk_ms, Re, diag] = compute_bulk_velocity( ...
        x_shifted, y_raw, U_flipped, D_int_px, D_ext_px, D_ext_mm, f_acq)
% COMPUTE_BULK_VELOCITY  Bulk velocity of a submerged jet from a PIV field.
%   Computes the volumetric-average velocity at the nozzle exit via
%   axisymmetric flow-rate integration: U_bulk = (Q/pi) / R^2,
%   where Q/pi = integral[ U(y) * |y - y_c| dy ].
%
% INPUTS:
%   x_shifted   PIV x-coordinates already mirrored and offset so that
%               x = 0 is at the nozzle exit plane [px]:
%               x_shifted = max(x_raw) - x_raw - nozzle_offset_px.
%   y_raw       PIV y-coordinates in the raw sensor frame [px].
%   U_flipped   PIV streamwise velocity with the sign flipped so that
%               positive values indicate flow in the jet direction.
%               U_flipped = -U_raw when the jet flows right-to-left.
%   D_int_px    Nozzle inner (flow) diameter [px]. Limits the integration
%               to the fluid region and excludes the nozzle wall.
%   D_ext_px    Nozzle outer diameter [px]. Used as the spatial length
%               scale for the px -> m conversion.
%   D_ext_mm    Nozzle outer diameter [mm]. Physical measurement from the
%               laboratory calibration.
%   f_acq       Acquisition frame rate [Hz]. Defines dt = 1/f_acq.
%
% OUTPUTS:
%   U_bulk      Bulk velocity [px/frame].
%   U_bulk_ms   Bulk velocity [m/s].
%   Re          Reynolds number: Re = U_bulk_ms * D_ext_m / nu_water,
%               where nu_water is the kinematic viscosity at 20 C.
%   diag        Diagnostic struct with the following fields:
%                 .y_center_px   Jet centroid position [px]
%                 .u_max         Peak velocity at the exit column [px/frame]
%                 .D_real_px     Detected jet flow diameter [px]
%                 .R_real_px     Detected jet flow radius [px]
%                 .x_salida      x-coordinate of the exit column [px]
%                 .px_to_m       Spatial scale factor [m/px]
%                 .dt            Frame interval [s/frame]

    if ~isequal(size(x_shifted), size(y_raw), size(U_flipped))
        error('compute_bulk_velocity: x_shifted, y_raw and U_flipped must have the same size.');
    end

    % --- Nozzle exit column: smallest positive x ---
    x_valid = x_shifted(x_shifted >= 0);
    if isempty(x_valid)
        error('compute_bulk_velocity: nozzle_offset_px is too large — no points with x >= 0.');
    end

    x_exit      = min(x_valid);
    mask_exit   = (x_shifted == x_exit);
    y_inlet     = y_raw(mask_exit);
    U_inlet     = U_flipped(mask_exit);

    % --- Jet centroid: velocity-weighted mean of the core (U > 50% of peak) ---
    U_sorted    = sort(U_inlet);
    idx_robust  = max(1, round(0.98 * length(U_sorted)));
    U_peak      = U_sorted(idx_robust);

    mask_core   = (U_inlet > 0.5 * U_peak);

    if sum(mask_core) < 2
        warning('compute_bulk_velocity: too few core points — centroid estimate unreliable.');
        y_center_px = mean(y_inlet);
    else
        % Velocity-weighted centroid (equivalent to the velocity centre of mass)
        y_center_px = sum(y_inlet(mask_core) .* U_inlet(mask_core)) / ...
                      sum(U_inlet(mask_core));
    end

    % --- Exit profile clipped to inner diameter ---
    mask_diam  = abs(y_raw - y_center_px) <= (D_int_px / 2);
    mask_valid = mask_exit & mask_diam;

    y_exit = y_raw(mask_valid);
    U_exit = U_flipped(mask_valid);

    if numel(y_exit) < 3
        error('compute_bulk_velocity: fewer than 3 points in the exit profile. Check D_int_px and nozzle_offset_px.');
    end

    % Sort spatially (required by trapz)
    [y_sorted, sort_idx] = sort(y_exit);
    u_sorted = U_exit(sort_idx);

    % --- Detect jet radius via 5% velocity threshold ---
    [u_max, idx_max] = max(u_sorted);
    y_center         = y_sorted(idx_max);   % peak location in the profile

    mask_fluid = u_sorted > (0.05 * u_max);
    y_fluid    = y_sorted(mask_fluid);
    u_fluid    = u_sorted(mask_fluid);

    if numel(y_fluid) < 3
        warning('compute_bulk_velocity: too few fluid points — result unreliable.');
        U_bulk = NaN;  U_bulk_ms = NaN;  Re = NaN;  diag = struct();
        return;
    end

    D_real_px = max(y_fluid) - min(y_fluid);
    R_real_px = D_real_px / 2;

    % --- Axisymmetric integration ---
    integrand = u_fluid .* abs(y_fluid - y_center);
    Q_over_pi = trapz(y_fluid, integrand);      % Q / pi  [px^3/frame]
    U_bulk    = Q_over_pi / (R_real_px^2);      % [px/frame]

    % --- Unit conversion ---
    px_to_m   = (D_ext_mm * 1e-3) / D_ext_px;   % [m/px]
    dt        = 1 / f_acq;                        % [s/frame]
    U_bulk_ms = U_bulk * (px_to_m / dt);          % [m/s]

    % --- Reynolds number (outer diameter, water at 20°C) ---
    nu_water  = 1.004e-6;           % [m^2/s]
    D_ext_m   = D_ext_mm * 1e-3;   % [m]
    Re        = U_bulk_ms * D_ext_m / nu_water;

    % --- Diagnostics ---
    diag.y_center_px = y_center_px;
    diag.u_max       = u_max;
    diag.D_real_px   = D_real_px;
    diag.R_real_px   = R_real_px;
    diag.x_salida    = x_exit;
    diag.px_to_m     = px_to_m;
    diag.dt          = dt;

    fprintf('\n--- compute_bulk_velocity ---\n');
    fprintf('  Exit column detected:            x = %.1f px\n',   x_exit);
    fprintf('  Jet centroid (Y):                    %.1f px\n',   y_center_px);
    fprintf('  Peak velocity (U_max):               %.3f px/frame\n', u_max);
    fprintf('  Detected flow diameter:              %.1f px\n',   D_real_px);
    fprintf('  Ratio U_max / U_bulk:                %.3f\n',      u_max / U_bulk);
    fprintf('  Spatial scale factor:                %.4f mm/px\n', px_to_m * 1e3);
    fprintf('  Frame interval (dt):                 %.2f ms\n',   dt * 1e3);
    fprintf('  U_bulk = %.3f px/frame  =  %.3f m/s\n',           U_bulk, U_bulk_ms);
    fprintf('  Reynolds number (Re):                %.0f\n',      Re);
    fprintf('-----------------------------\n\n');

end