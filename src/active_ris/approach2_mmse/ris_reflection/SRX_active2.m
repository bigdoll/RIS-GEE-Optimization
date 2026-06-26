function srx_active = SRX_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
% SRX_ACTIVE2  MMSE sum-rate in X as F1(X) - F2(X).
%
%   function srx_active = SRX_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (50)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    G1 = G1func_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW);
    G2 = G2func_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW);
   
    srx_active = G1 - G2;
    
end
