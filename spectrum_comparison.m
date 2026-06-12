% SPECTRUM_COMPARISON
% Overlays Welch PSD of PaIRS and a custom algorithm at two measurement
% points (Shear Layer and Jet Core) for each particle density dataset.
% Requires data_point_selection.m run twice: once with source='PaIRS'
% and once with source='<algorithm>'.

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
factor_N           = 8;
algo_name          = 'SDC Predictor';
algo_tag           = matlab.lang.makeValidName(algo_name);
ref_tag            = 'PaIRS';
D_ext_m            = 0.030;                    % [m]
U_bulk_per_density = [0.148, 0.150];           % [m/s]
color_ref          = [0.10 0.10 0.45];         % PaIRS (dark blue)
color_algo         = [0.75 0.10 0.10];         % algorithm (dark red)


%% 2. MAIN LOOP
for d = 1:length(densities)
    dens   = densities{d};
    U_bulk = U_bulk_per_density(d);

    f_SL_ref    = sprintf('data_spectrum_ShearLayer_%s_%s.mat', dens, ref_tag);
    f_Core_ref  = sprintf('data_spectrum_Core_%s_%s.mat',       dens, ref_tag);
    f_SL_algo   = sprintf('data_spectrum_ShearLayer_%s_%s.mat', dens, algo_tag);
    f_Core_algo = sprintf('data_spectrum_Core_%s_%s.mat',       dens, algo_tag);

    if ~all(cellfun(@(f) exist(f,'file'), {f_SL_ref, f_Core_ref, f_SL_algo, f_Core_algo}))
        fprintf('Skipping %s: missing .mat files.\n', dens);
        continue;
    end

    SL_ref    = load(f_SL_ref,    'v_time_SL',   'fs', 'x_SL_target_norm', 'y_SL_target_norm');
    Core_ref  = load(f_Core_ref,  'u_time_Core', 'fs', 'x_Core_target_norm');
    SL_algo   = load(f_SL_algo,   'v_time_SL',   'fs');
    Core_algo = load(f_Core_algo, 'u_time_Core', 'fs');
    fs = SL_ref.fs;

    % Truncate to shortest length for identical frequency resolution
    N = min([length(SL_ref.v_time_SL), length(SL_algo.v_time_SL), ...
             length(Core_ref.u_time_Core), length(Core_algo.u_time_Core)]);
    SL_ref.v_time_SL      = SL_ref.v_time_SL(1:N);
    SL_algo.v_time_SL     = SL_algo.v_time_SL(1:N);
    Core_ref.u_time_Core  = Core_ref.u_time_Core(1:N);
    Core_algo.u_time_Core = Core_algo.u_time_Core(1:N);

    % --- Welch PSD ---
    window_size = floor(N / factor_N);
    overlap     = floor(window_size / 2);   % 50% overlap

    [P_SL_ref,    f] = pwelch(SL_ref.v_time_SL    - mean(SL_ref.v_time_SL),       window_size, overlap, window_size, fs);
    [P_SL_algo,   ~] = pwelch(SL_algo.v_time_SL   - mean(SL_algo.v_time_SL),      window_size, overlap, window_size, fs);
    [P_Core_ref,  ~] = pwelch(Core_ref.u_time_Core  - mean(Core_ref.u_time_Core),  window_size, overlap, window_size, fs);
    [P_Core_algo, ~] = pwelch(Core_algo.u_time_Core - mean(Core_algo.u_time_Core), window_size, overlap, window_size, fs);

    % --- Strouhal axis ---
    St = f * D_ext_m / U_bulk;

    % --- Kolmogorov -5/3 reference: anchored to PaIRS at 15% of St range ---
    idx_ref  = floor(length(St)*0.15) : floor(length(St)*0.5);
    St_ref   = St(idx_ref);
    P_k_SL   = P_SL_ref(idx_ref(1))   * (St_ref(1)^(5/3)) * St_ref.^(-5/3);
    P_k_Core = P_Core_ref(idx_ref(1)) * (St_ref(1)^(5/3)) * St_ref.^(-5/3);

    % --- Figure A: Shear Layer ---
    fig_SL = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 550]);
    loglog(St, P_SL_ref,  'Color', color_ref,  'LineWidth', 1.6, 'DisplayName', "PaIRS (v')"); hold on;
    loglog(St, P_SL_algo, 'Color', color_algo, 'LineWidth', 1.6, 'DisplayName', sprintf("%s (v')", algo_name));
    loglog(St_ref, P_k_SL, 'k--', 'LineWidth', 1.8, 'DisplayName', 'Slope -5/3 (Kolmogorov)');
    setup_plot(FONT_SIZE); legend('Location', 'southwest', 'FontSize', LEGEND_SIZE);
    save_fig(fig_SL, output_dir, dens, ['ShearLayer_' algo_tag], factor_N);

    % --- Figure B: Jet Core ---
    fig_Core = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 550]);
    loglog(St, P_Core_ref,  'Color', color_ref,  'LineWidth', 1.6, 'DisplayName', "PaIRS (u')"); hold on;
    loglog(St, P_Core_algo, 'Color', color_algo, 'LineWidth', 1.6, 'DisplayName', sprintf("%s (u')", algo_name));
    loglog(St_ref, P_k_Core, 'k--', 'LineWidth', 1.8, 'DisplayName', 'Slope -5/3 (Kolmogorov)');
    setup_plot(FONT_SIZE); legend('Location', 'southwest', 'FontSize', LEGEND_SIZE);
    save_fig(fig_Core, output_dir, dens, ['Core_' algo_tag], factor_N);

    fprintf('Done: %s vs PaIRS | %s | U_bulk=%.3f m/s\n', algo_name, dens, U_bulk);
end

disp('--- Spectral comparison complete ---');


%% HELPER FUNCTIONS
function setup_plot(font_sz)
    grid on; grid minor;
    set(gca, 'TickLabelInterpreter', 'tex', 'FontName', 'Times New Roman', 'FontSize', font_sz);
    xlabel('St = f \cdot D / U_{bulk}', 'Interpreter', 'tex', 'FontSize', font_sz);
    ylabel('PSD \Phi  (m^{2} s^{-2} Hz^{-1})', 'Interpreter', 'tex', 'FontSize', font_sz);
end

function save_fig(fig_obj, output_dir, dens, tag, factor)
    exportgraphics(fig_obj, fullfile(output_dir, ...
        sprintf('SpectrumComp_%s_%s_WinN%d.png', dens, tag, factor)), 'Resolution', 300);
    close(fig_obj);
end