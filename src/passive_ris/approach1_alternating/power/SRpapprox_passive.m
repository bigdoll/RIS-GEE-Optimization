function srp_approx = SRpapprox_passive(p, p_bar, G, H, gamma, sigma_sq, K, NR, BW)
% SRPAPPROX_PASSIVE  Concave surrogate of the sum-rate in p:  g1 - g2 - grad g2'(p-p_bar).
%
%   function srp_approx = SRpapprox_passive(p, p_bar, G, H, gamma, sigma_sq, K, NR, BW)
%
% Paper reference: Eq. (41)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    % Linear MMSE update
    C = LMMSE_receiver_passive(p_bar, G, H, gamma, sigma_sq, K, NR);

    g1 = func_g1(p, C, G, H, gamma, sigma_sq, K, BW);
    g2 = func_g2(p_bar, C, G, H, gamma, sigma_sq, K, BW);
    grad_g2 = gradient_g2(p_bar, C, G, H, gamma, sigma_sq, K, BW);
    
    srp_approx = g1 - g2 - grad_g2'*(p - p_bar);
    
end
