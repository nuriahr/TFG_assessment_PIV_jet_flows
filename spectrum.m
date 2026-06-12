% SPECTRUM
% Welch PSD of velocity fluctuations at the Shear Layer (v') and Jet Core (u').
% Generates three figures per density: comparative, Shear Layer only, Core only.
% Requires data_point_selection.m to be run first (source = 'PaIRS').

clear; clc; close all;

%% 1. CONFIGURATION
output_dir = 'path/to/output/Spectrum';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

FONT_SIZE   = 26;
LEGEND_SIZE = 18;

set(groot, 'defaultAxesFontName',             'Times New Roman');
set(groot, 'defaultTextFontName',             'Times New Roman');
set(groot, 'defaultLegendFontName',           'Times New Roman');
set(groot, 'defaultAxesFontSize',             FONT_SIZE);
set(groot, 'defaultTextFontSize',             FONT_SIZE);
set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
set(groot, 'defaultTextInterpreter',          'tex');
set(groot, 'defaultLegendInterpreter',        'tex');

densities          = {'medium', 'high'};
factor_N           = 16;          % Welch window = N / factor_N
D_ext_m            = 0.030;       % [m]
U_bulk_per_density = [0.148, 0.150];   % [m/s], same order as densities


%% 2. MAIN LOOP
for d = 1:length(densities)
    dens   = densities{d};
    U_bulk = U_bulk_per_density(d);

    file_SL   = sprintf('data_spectrum_ShearLayer_%s_PaIRS.mat', dens);
    file_Core = sprintf('data_spectrum_Core_%s_PaIRS.mat', dens);
    if ~exist(file_SL, 'file') || ~exist(file_Core, 'file')
        fprintf('Skipping %s: input .mat files not found.\n', dens);
        continue;
    end

    data_SL   = load(file_SL,   'v_time_SL',   'fs', 'x_SL_target_norm',   'y_SL_target_norm');
    data_Core = load(file_Core, 'u_time_Core', 'fs', 'x_Core_target_norm', 'y_Core_target_norm');
    fs = data_SL.fs;

    % --- Welch PSD ---
    N           = length(data_SL.v_time_SL);
    window_size = floor(N / factor_N);
    overlap     = floor(window_size / 2);   % 50% overlap

    df = fs / window_size;
    fprintf('%s | N=%d | df=%.2f Hz | St=[%.3f, %.2f]\n', upper(dens), N, df, ...
            df * D_ext_m / U_bulk, (fs/2) * D_ext_m / U_bulk);

    [P_SL,   f] = pwelch(data_SL.v_time_SL   - mean(data_SL.v_time_SL),   window_size, overlap, window_size, fs);
    [P_Core, ~] = pwelch(data_Core.u_time_Core - mean(data_Core.u_time_Core), window_size, overlap, window_size, fs);

    % --- Strouhal axis ---
    St = f * D_ext_m / U_bulk;

    % --- Kolmogorov -5/3 reference: anchored at 15% of St range ---
    idx_ref = floor(length(St)*0.15) : floor(length(St)*0.5);
    St_ref  = St(idx_ref);
    P_kolmogorov_SL   = P_SL(idx_ref(1))   * (St_ref(1)^(5/3)) * St_ref.^(-5/3);
    P_kolmogorov_Core = P_Core(idx_ref(1)) * (St_ref(1)^(5/3)) * St_ref.^(-5/3);

    lbl_SL   = sprintf('Shear Layer (x/D=%.1f, y/D=%.1f, v'')', data_SL.x_SL_target_norm,    data_SL.y_SL_target_norm);
    lbl_Core = sprintf('Jet Core (x/D=%.1f, y/D=%.1f, u'')',    data_Core.x_Core_target_norm, data_Core.y_Core_target_norm);

    % --- Figure A: Comparative ---
    fig_comp = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 550]);
    loglog(St, P_SL, 'r', 'LineWidth', 1.5, 'DisplayName', lbl_SL); hold on;
    loglog(St, P_Core, 'b', 'LineWidth', 1.5, 'DisplayName', lbl_Core);
    loglog(St_ref, P_kolmogorov_SL, 'k--', 'LineWidth', 1.8, 'DisplayName', 'Slope -5/3 (Kolmogorov)');
    setup_plot(FONT_SIZE); legend('Location', 'southwest', 'FontSize', LEGEND_SIZE);
    save_fig(fig_comp, output_dir, dens, 'Comparison', factor_N);

    % --- Figure B: Shear Layer ---
    fig_SL = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 550]);
    loglog(St, P_SL, 'r', 'LineWidth', 1.5, 'DisplayName', lbl_SL); hold on;
    loglog(St_ref, P_kolmogorov_SL, 'k--', 'LineWidth', 1.8, 'DisplayName', 'Slope -5/3 (Kolmogorov)');
    setup_plot(FONT_SIZE); legend('Location', 'southwest', 'FontSize', LEGEND_SIZE);
    save_fig(fig_SL, output_dir, dens, 'ShearLayer', factor_N);

    % --- Figure C: Jet Core ---
    fig_core = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 550]);
    loglog(St, P_Core, 'b', 'LineWidth', 1.5, 'DisplayName', lbl_Core); hold on;
    loglog(St_ref, P_kolmogorov_Core, 'k--', 'LineWidth', 1.8, 'DisplayName', 'Slope -5/3 (Kolmogorov)');
    setup_plot(FONT_SIZE); legend('Location', 'southwest', 'FontSize', LEGEND_SIZE);
    save_fig(fig_core, output_dir, dens, 'Core', factor_N);

    fprintf('Done: %s | N/%d | U_bulk=%.3f m/s\n', upper(dens), factor_N, U_bulk);
end

disp('--- Spectral analysis complete ---');


%% HELPER FUNCTIONS
function setup_plot(font_sz)
    grid on; grid minor;
    set(gca, 'TickLabelInterpreter', 'tex', 'FontName', 'Times New Roman', 'FontSize', font_sz);
    xlabel('St = f \cdot D / U_{bulk}', 'Interpreter', 'tex', 'FontSize', font_sz);
    ylabel('PSD \Phi  (m^{2} s^{-2} Hz^{-1})', 'Interpreter', 'tex', 'FontSize', font_sz);
end

function save_fig(fig_obj, output_dir, dens, tag, factor)
    exportgraphics(fig_obj, fullfile(output_dir, ...
        sprintf('Spectrum_%s_%s_WinN%d.png', dens, tag, factor)), 'Resolution', 300);
    close(fig_obj);
end