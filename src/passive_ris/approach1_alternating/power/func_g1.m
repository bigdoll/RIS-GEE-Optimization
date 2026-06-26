function g1 = func_g1(p, C, G, H, gamma, sigma_sq, K, BW)
% FUNC_G1  Concave term g1(p) of the sum-rate decomposition.
%
%   function g1 = func_g1(p, C, G, H, gamma, sigma_sq, K, BW)
%
% Paper reference: Eq. (40)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    g1 = 0; % initialization
  
    for k=1:K
        ck = C(:,k); % Linear MMSE received filter for user k 
        pakk_sum = 0;
        
        for m=1:K
            hm = H(:,m); % channel for user m b/w UE-RIS
            Hm = diag(hm); % diagonal matrix for hm
            Am = G * Hm; % Cascaded channel matrix for user m
            akm = abs(ck'*Am*gamma)^2; 
            pm = p(m)/sigma_sq;
            pakk_sum = pakk_sum + pm*akm; % sum of all the signal power 
        end
        dk = (norm(ck))^2; %sigma_sq*(norm(ck))^2;
        g1 = g1 + BW*log2(sigma_sq) + BW*log2(dk + pakk_sum);
        
    end
   
end
