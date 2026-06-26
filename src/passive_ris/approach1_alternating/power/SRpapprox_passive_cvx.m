function srpapprox_cvx = SRpapprox_passive_cvx(p, p_bar, G, H, gamma, sigma_sq, K, NR, BW)
% SRPAPPROX_PASSIVE_CVX  CVX-expression form of the power surrogate.
%
%   function srpapprox_cvx = SRpapprox_passive_cvx(p, p_bar, G, H, gamma, sigma_sq, K, NR, BW)
%
% Paper reference: Eq. (41)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    % Linear MMSE update
    C = LMMSE_receiver_passive(p_bar, G, H, gamma, sigma_sq, K, NR);

    g1_cvx = func_g1_cvx(p, C, G, H, gamma, sigma_sq, K, BW);
    grad_g2_cvx = gradient_g2_cvx(p_bar, C, G, H, gamma, sigma_sq, K, BW);
    
    srpapprox_cvx = g1_cvx - grad_g2_cvx' * p; % sequential approximation for sum_rate
    
end
