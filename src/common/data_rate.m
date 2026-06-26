function [user_rate, sum_rate] = data_rate(sinr, K, BW)
% DATA_RATE  Per-user and sum rate  R_k = B*log2(1+SINR_k).
%
%   function [user_rate, sum_rate] = data_rate(sinr, K, BW)
%
% Paper reference: Eq. (23)/(123)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    % Initialization
    user_rate = zeros(K,1);
    sum_rate = 0;
    
    for k=1:K
        user_rate(k) = BW * log2(1 + sinr(k));
        sum_rate = sum_rate + user_rate(k);
    end
    
end
