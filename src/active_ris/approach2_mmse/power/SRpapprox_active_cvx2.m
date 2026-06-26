function srpapprox_cvx = SRpapprox_active_cvx2(p, p_bar, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
% SRPAPPROX_ACTIVE_CVX2  CVX-expression form of the power surrogate (MMSE branch).
%
%   function srpapprox_cvx = SRpapprox_active_cvx2(p, p_bar, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (55)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    g1a_cvx = funcg1_active_cvx2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW);
    grad_g2a = gradg2_active2(p_bar, G, H, X, sigma_sq, sigma_sqris, K, NR, BW);

    srpapprox_cvx = g1a_cvx - grad_g2a' * p; % sequential approximation for sum_rate
    
end
