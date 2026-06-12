function [RMS_error, RMS_core, E_map] = compute_piv_error(X_ref, Y_ref, U_ref, V_ref, ...
                                                          X_eval, Y_eval, U_eval, V_eval, ...
                                                          core_mask)
% COMPUTE_PIV_ERROR  Velocity magnitude error metrics vs. ground truth.
%   Interpolates the evaluated field onto the reference grid and returns
%   the global RMS error, the core-region RMS, and the spatial error map.
%
% INPUTS:
%   X_ref, Y_ref, U_ref, V_ref      Reference field (PaIRS / ground truth).
%   X_eval, Y_eval, U_eval, V_eval  Evaluated field (algorithm output).
%   core_mask  (optional) Logical mask, same size as the reference grid.
%              true = cells belonging to the central jet-core region.
%              If omitted, RMS_core is returned as NaN.
%
% OUTPUTS:
%   RMS_error  Global RMS of the magnitude error over the common region [px/frame].
%   RMS_core   RMS of the magnitude error restricted to core_mask [px/frame].
%              Returns NaN if core_mask is not provided or has no valid cells.
%   E_map      Spatial map of the magnitude error [px/frame].
%              NaN outside the common validity region.

    if nargin < 9
        core_mask = [];
    end

    % --- Interpolate evaluated field onto reference grid ---
    % NaN entries (unmeasured windows) are zeroed before interpolation;
    % the common mask below excludes them from all metrics.
    U_eval_clean = U_eval;  U_eval_clean(isnan(U_eval)) = 0;
    V_eval_clean = V_eval;  V_eval_clean(isnan(V_eval)) = 0;

    U_eval_interp = interp2(X_eval, Y_eval, U_eval_clean, X_ref, Y_ref, 'spline', 0);
    V_eval_interp = interp2(X_eval, Y_eval, V_eval_clean, X_ref, Y_ref, 'spline', 0);

    % --- Common validity mask: non-zero, non-NaN in both fields ---
    mask_ref    = (U_ref ~= 0) & (V_ref ~= 0) & ~isnan(U_ref);
    mask_eval   = (U_eval_interp ~= 0) & (V_eval_interp ~= 0) & ~isnan(U_eval_interp);
    mask_common = mask_ref & mask_eval;

    % --- Magnitude error and global RMS ---
    Mag_ref  = sqrt(U_ref.^2  + V_ref.^2);
    Mag_eval = sqrt(U_eval_interp.^2 + V_eval_interp.^2);

    E_mag             = abs(Mag_eval - Mag_ref);
    E_map             = E_mag;
    E_map(~mask_common) = NaN;   % undefined outside the common region

    if sum(mask_common(:)) > 0
        RMS_error = sqrt(mean(E_mag(mask_common).^2));
    else
        RMS_error = NaN;
    end

     % --- Jet core RMS ---
    RMS_core = NaN;
    if ~isempty(core_mask)
        mask_core_valid = mask_common & core_mask;
        if sum(mask_core_valid(:)) > 0
            RMS_core = sqrt(mean(E_mag(mask_core_valid).^2));
        end
    end
end