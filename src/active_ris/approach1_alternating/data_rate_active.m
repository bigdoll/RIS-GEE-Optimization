function [ur_a, sr_a] = data_rate_active(sinr_a, K, BW)
% DATA_RATE_ACTIVE  Per-user and sum rate for the active branch (bandwidth in Mbit/s).
%
%   function [ur_a, sr_a] = data_rate_active(sinr_a, K, BW)
%
% Paper reference: Eq. (23)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    % Initialization
    ur_a = zeros(K,1);
    sr_a = 0;
    
    for k=1:K
        ur_a(k) = BW * log2(1 + sinr_a(k));
        sr_a = sr_a + ur_a(k);
    end
    
end
