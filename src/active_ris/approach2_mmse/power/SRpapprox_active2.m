function srp_approx = SRpapprox_active2(p, p_bar, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
% SRPAPPROX_ACTIVE2  Concave surrogate of the sum-rate in p (MMSE-embedded branch).
%
%   function srp_approx = SRpapprox_active2(p, p_bar, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (55)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    g1 = funcg1_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW);
    g2 = funcg2_active2(p_bar, G, H, X, sigma_sq, sigma_sqris, K, NR, BW);
    grad_g2 = gradg2_active2(p_bar, G, H, X, sigma_sq, sigma_sqris, K, NR, BW);
    
    srp_approx = g1 - g2 - grad_g2'*(p - p_bar);
    
end
