function srp_active = SRp_active(p, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW)
% SRP_ACTIVE  MMSE sum-rate written as g1(p) - g2(p).
%
%   function srp_active = SRp_active(p, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (40)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    % Linear MMSE update
    C = LMMSE_receiver_active(p, G, H, gamma, sigma_sq, sigma_sqris, K, NR);

    g1 = funcg1_active(p, C, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW);
    g2 = funcg2_active(p, C, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW);
   
    srp_active = g1 - g2;
    
end
