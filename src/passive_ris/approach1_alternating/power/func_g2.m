function g2 = func_g2(p, C, G, H, gamma, sigma_sq, K, BW)
% FUNC_G2  Convex term g2(p) of the sum-rate decomposition.
%
%   function g2 = func_g2(p, C, G, H, gamma, sigma_sq, K, BW)
%
% Paper reference: Eq. (40)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    g2 = 0; % initialization
    
    for k=1:K
        ck = C(:,k); % Linear MMSE received filter for user k 
        pakm_sum = 0;
        
        for m=1:K
            if m ~= k
                hm = H(:,m); % channel for user m b/w UE-RIS
                Hm = diag(hm); % diagonal matrix for hm
                Am = G * Hm; % Cascaded channel matrix for user m
                akm = abs(ck'*Am*gamma)^2;
                pm = p(m)/sigma_sq;
                pakm_sum = pakm_sum + pm*akm; % sum of all the signal power 
            end
        end
        
        dk = (norm(ck))^2;
        g2 = g2 + BW*log2(sigma_sq) + BW*log2(dk + pakm_sum);
        
    end
    
end
