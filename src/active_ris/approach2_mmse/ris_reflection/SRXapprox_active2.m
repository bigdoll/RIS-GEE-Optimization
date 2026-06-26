function srX_approx = SRXapprox_active2(p, G, H, X, X_bar, sigma_sq, sigma_sqris, K, NR, BW)
% SRXAPPROX_ACTIVE2  Concave surrogate of the sum-rate in X (F2 linearized).
%
%   function srX_approx = SRXapprox_active2(p, G, H, X, X_bar, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (52)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    G1 = G1func_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW);
    G2 = G2func_active2(p, G, H, X_bar, sigma_sq, sigma_sqris, K, NR, BW);
    grad_G2 = G2grad_active2(p, G, H, X_bar, sigma_sq, sigma_sqris, K, NR, BW);
    
    srX_approx = G1 - G2 - real(trace(grad_G2'*(X - X_bar)));
    
end
