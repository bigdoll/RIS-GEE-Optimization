function G1a_cvx = G1func_active_cvx2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
% G1FUNC_ACTIVE_CVX2  CVX-expression form of F1(X).
%
%   function G1a_cvx = G1func_active_cvx2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (50)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    G1a_cvx = 0; % initialization
    N =  length(X);
    diag_X = X.*eye(N);
    sigma_norm = sigma_sqris / sigma_sq;
    W_norm = eye(NR) + sigma_norm*(G*diag_X*G');
%     W = sigma_sq*eye(NR) + sigma_sqris*(G*diag_X*G');

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

        G1a_cvx = G1a_cvx + BW*log_det(W_norm + pakk_sum)/log(2);
        
    end
   
end
