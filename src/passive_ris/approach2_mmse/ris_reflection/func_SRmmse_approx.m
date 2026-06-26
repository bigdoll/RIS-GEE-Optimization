function SRapprox_mmse = func_SRmmse_approx(p, G, H, X, X_bar, sigma_sq, K, NR, BW)
% FUNC_SRMMSE_APPROX  Concave surrogate of the sum-rate in X (F2 linearized).
%
%   function SRapprox_mmse = func_SRmmse_approx(p, G, H, X, X_bar, sigma_sq, K, NR, BW)
%
% Paper reference: Eq. (52)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    [G1, ~] = func_G(p, G, H, X, sigma_sq, K, NR, BW);
    [~, G2] = func_G(p, G, H, X_bar, sigma_sq, K, NR, BW);
    grad_G2 = funcG2_grad(p, G, H, X_bar, sigma_sq, K, NR, BW);
    G2_approx = G2 + real(trace(grad_G2'*(X-X_bar)));
    
    SRapprox_mmse = G1 - G2_approx;
end
