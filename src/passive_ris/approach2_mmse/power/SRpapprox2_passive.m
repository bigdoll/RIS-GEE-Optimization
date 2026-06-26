function srp2_approx = SRpapprox2_passive(p, p_bar, G, H, X, sigma_sq, K, NR, BW, option)
% SRPAPPROX2_PASSIVE  Concave surrogate of the sum-rate in p (MMSE-embedded branch).
%
%   function srp2_approx = SRpapprox2_passive(p, p_bar, G, H, X, sigma_sq, K, NR, BW, option)
%
% Paper reference: Eq. (55)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
   
    f1 = func2_g1(p, G, H, X, sigma_sq, K, NR, BW, option);
    f2 = func2_g2(p_bar, G, H, X, sigma_sq, K, NR, BW, option);
    grad_f2 = func2_gradient(p_bar, G, H, X, sigma_sq, K, NR, BW);
    
    srp2_approx = f1 - f2 - grad_f2'*(p - p_bar);
    
end
