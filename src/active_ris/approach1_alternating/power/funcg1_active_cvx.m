function g1a_cvx = funcg1_active_cvx(p, C, G, H, gamma, sigma_sq, sigma_sqris, K,  NR, BW)
% FUNCG1_ACTIVE_CVX  CVX-expression form of g1(p).
%
%   function g1a_cvx = funcg1_active_cvx(p, C, G, H, gamma, sigma_sq, sigma_sqris, K,  NR, BW)
%
% Paper reference: Eq. (40)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    g1a_cvx = 0; % initialization
    Gamma = diag(gamma);
    W = sigma_sq*eye(NR) + sigma_sqris*G*(Gamma*Gamma')*G';

    for k=1:K
        ck = C(:,k); % Linear MMSE received filter for user k 
        pakk_sum = 0;
        dk = real(ck'*W*ck);

        for m=1:K
            hm = H(:,m); % channel for user m b/w UE-RIS
            Hm = diag(hm); % diagonal matrix for hm
            Am = G * Hm; % Cascaded channel matrix for user m
            akm = abs(ck'*Am*gamma)^2;
            pm = p(m)/dk;
            pakk_sum = pakk_sum + pm*akm; % sum of all the signal power 
        end 

        g1a_cvx = g1a_cvx + BW*log(1 + pakk_sum)/log(2); % BW*log(sigma_sq*norm(ck)^2)/log(2)
        
    end
   
end
