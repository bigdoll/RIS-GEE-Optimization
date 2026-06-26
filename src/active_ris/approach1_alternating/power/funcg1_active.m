function g1a = funcg1_active(p, C, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW)
% FUNCG1_ACTIVE  Concave term g1(p) of the sum-rate decomposition.
%
%   function g1a = funcg1_active(p, C, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (40)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    g1a = 0; % initialization
    Gamma = diag(gamma);
    W = sigma_sq*eye(NR) + sigma_sqris*G*(Gamma*Gamma')*G';

    for k=1:K
        ck = C(:,k); % Linear MMSE received filter for user k 
        pakk_sum = 0;
      
        for m=1:K
            hm = H(:,m); % channel for user m b/w UE-RIS
            Hm = diag(hm); % diagonal matrix for hm
            Am = G * Hm; % Cascaded channel matrix for user m
            akm = abs(ck'*Am*gamma)^2;
            pm = p(m);
            pakk_sum = pakk_sum + pm*akm; % sum of all the signal power 
        end 
        dk = real(ck'*W*ck);
        g1a = g1a + BW*log(dk + pakk_sum)/log(2); % BW*log(sigma_sq*norm(ck)^2)/log(2)
        
    end
   
end
