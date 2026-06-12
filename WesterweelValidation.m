function [U_out, V_out] = WesterweelValidation(U_in, V_in, Thr, eps_val)
% WESTERWEELVALIDATION  Universal outlier detection (Westerweel & Scarano 2005).
%   Detects outliers using the normalised median residual test applied
%   jointly to U and V, and replaces them with the local median.
%
% INPUTS:
%   U_in, V_in     Input velocity fields
%   Thr            Detection threshold. Recommended range: 2.0–3.0.
%   eps_val        Noise floor (epsilon). Typically 0.1.
%
% OUTPUTS:
%   U_out, V_out   Output velocity fields

    [J, I] = size(U_in);
    
    % Initialize outputs
    U_out = U_in;
    V_out = V_in;
    
    Median_U = zeros(J, I);
    Median_V = zeros(J, I);
    NormFluct_U = zeros(J, I);
    NormFluct_V = zeros(J, I);
    
    b = 1; % neighbourhood half-width (3x3 window)
    
    % --- Step 1: normalised fluctuation for each interior vector ---
    for i = (1+b):(I-b)  
        for j = (1+b):(J-b)  
            
            % --- U component ---
            Neigh_U = U_in(j-b:j+b, i-b:i+b); 
            NeighCol_U = Neigh_U(:);
            NeighCol2_U = NeighCol_U([1:4, 6:9]); % exclude centre

            med_u = median(NeighCol2_U);
            Median_U(j,i) = med_u;
            
            % Residual and fluctuation
            res_u = abs(NeighCol2_U - med_u);
            median_res_u = median(res_u);
            fluct_u = abs(U_in(j,i) - med_u);
            
            % Normalised fluctuation
            NormFluct_U(j,i) = fluct_u / (median_res_u + eps_val);
            
            % --- V component (same process) ---
            Neigh_V = V_in(j-b:j+b, i-b:i+b);
            NeighCol_V = Neigh_V(:);
            NeighCol2_V = NeighCol_V([1:4, 6:9]);
            
            med_v = median(NeighCol2_V);
            Median_V(j,i) = med_v;
            
            res_v = abs(NeighCol2_V - med_v);
            median_res_v = median(res_v);
            fluct_v = abs(V_in(j,i) - med_v);
            
            NormFluct_V(j,i) = fluct_v / (median_res_v + eps_val);
            
        end
    end
    
    % --- Step 2: combined vector residual and outlier mask --- 
    Combined_Res = sqrt(NormFluct_U.^2 + NormFluct_V.^2);
    Outlier_Mask = Combined_Res > Thr;
    
    % --- Step 3: replace outliers with the local median ---
    % Index of outlier vectors
    idx = find(Outlier_Mask);

    % Apply corrections
    U_out(idx) = Median_U(idx);
    V_out(idx) = Median_V(idx);
    
end