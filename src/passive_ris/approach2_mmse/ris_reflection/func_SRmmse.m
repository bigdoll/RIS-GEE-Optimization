function SR_mmse = func_SRmmse(p, G, H, X, sigma_sq, K, NR, BW)
% FUNC_SRMMSE  MMSE sum-rate in X as F1(X) - F2(X).
%
%   function SR_mmse = func_SRmmse(p, G, H, X, sigma_sq, K, NR, BW)
%
% Paper reference: Eq. (50)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    [G1, G2] = func_G(p, G, H, X, sigma_sq, K, NR, BW);
    
    SR_mmse = G1 - G2;
end
