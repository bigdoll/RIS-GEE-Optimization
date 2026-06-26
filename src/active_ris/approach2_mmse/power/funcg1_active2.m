function g1a = funcg1_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
% FUNCG1_ACTIVE2  Concave term g1(p) of the MMSE-embedded sum-rate.
%
%   function g1a = funcg1_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (55)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    g1a = 0; % initialization
    gamma_vec = diag(X, 0);
    Gamma_herm = diag(gamma_vec); % Gamma_herm = Gamma*Gamma'
    W = sigma_sq*eye(NR) + sigma_sqris*G*Gamma_herm*G';

    for k=1:K
        pakk_sum = 0;
      
        for m=1:K
            hm = H(:,m); % channel for user m b/w UE-RIS
            Hm = diag(hm); % diagonal matrix for hm
            Am = G * Hm; % Cascaded channel matrix for user m
            akm = Am*X*Am';
            pm = p(m);
            pakk_sum = pakk_sum + pm*akm; % sum of all the signal power 
        end

        g1a = g1a + BW*log(real(det(W + pakk_sum)))/log(2);
        
    end
   
end
