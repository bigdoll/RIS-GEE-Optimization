function g2a = funcg2_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
% FUNCG2_ACTIVE2  Convex term g2(p) of the MMSE-embedded sum-rate.
%
%   function g2a = funcg2_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (55)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    g2a = 0; % initialization
    gamma_vec = diag(X, 0);
    Gamma_herm = diag(gamma_vec); % Gamma_herm = Gamma*Gamma'
    W = sigma_sq*eye(NR) + sigma_sqris*G*Gamma_herm*G';
    
    for k=1:K
        pakm_sum = 0;
        
        for m=1:K
            if m ~= k
                hm = H(:,m); % channel for user m b/w UE-RIS
                Hm = diag(hm); % diagonal matrix for hm
                Am = G * Hm; % Cascaded channel matrix for user m
                akm = Am*X*Am';
                pm = p(m);
                pakm_sum = pakm_sum + pm*akm; % sum of all the signal power
            end
        end

        g2a = g2a + BW*log(real(det(W + pakm_sum)))/log(2); 
        
    end
   
end
