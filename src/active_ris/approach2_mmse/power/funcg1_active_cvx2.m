function g1a_cvx2 = funcg1_active_cvx2(p, G, H, X, sigma_sq, sigma_sqris, K,  NR, BW)
% FUNCG1_ACTIVE_CVX2  CVX-expression form of g1(p) (MMSE branch).
%
%   function g1a_cvx2 = funcg1_active_cvx2(p, G, H, X, sigma_sq, sigma_sqris, K,  NR, BW)
%
% Paper reference: Eq. (55)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    g1a_cvx2 = 0; % initialization
    gamma_vec = diag(X, 0);
    Gamma_herm = diag(gamma_vec); % Gamma_herm = Gamma*Gamma'
%     W = sigma_sq*eye(NR) + sigma_sqris*G*(Gamma_herm)*G';
    sigma_norm = sigma_sqris/sigma_sq;

    W_norm = eye(NR) + sigma_norm*(G*Gamma_herm*G');

    for k=1:K
        pakk_sum = 0;

        for m=1:K
            hm = H(:,m); % channel for user m b/w UE-RIS
            Hm = diag(hm); % diagonal matrix for hm
            Am = G * Hm; % Cascaded channel matrix for user m
            akm = Am*X*Am';
            pm = p(m)/sigma_sq;
            pakk_sum = pakk_sum + pm*akm; % sum of all the signal power 
        end 

        g1a_cvx2 = g1a_cvx2 + BW*log_det(W_norm + pakk_sum)/log(2); % BW*log(sigma_sq*norm(ck)^2)/log(2)
        
    end
   
end
