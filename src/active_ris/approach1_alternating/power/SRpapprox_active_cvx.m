function srpapprox_cvx = SRpapprox_active_cvx(p, p_bar, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW)
% SRPAPPROX_ACTIVE_CVX  CVX-expression form of the power surrogate.
%
%   function srpapprox_cvx = SRpapprox_active_cvx(p, p_bar, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (41)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    % Linear MMSE update
    C = LMMSE_receiver_active(p_bar, G, H, gamma, sigma_sq, sigma_sqris, K, NR);

    g1a_cvx = funcg1_active_cvx(p, C, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW);
    grad_g2a = gradg2_active(p_bar, C, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW);
    
    srpapprox_cvx = g1a_cvx - grad_g2a' * p; % sequential approximation for sum_rate
    
end
