function srp_approx = SRpapprox_active(p, p_bar, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW)
% SRPAPPROX_ACTIVE  Concave surrogate of the sum-rate in p:  g1 - g2 - grad g2'(p-p_bar).
%
%   function srp_approx = SRpapprox_active(p, p_bar, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (41)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    % Linear MMSE update
    C = LMMSE_receiver_active(p_bar, G, H, gamma, sigma_sq, sigma_sqris, K, NR);

    g1 = funcg1_active(p, C, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW);
    g2 = funcg2_active(p_bar, C, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW);
    grad_g2 = gradg2_active(p_bar, C, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW);
    
    srp_approx = g1 - g2 - grad_g2'*(p - p_bar);
    
end
