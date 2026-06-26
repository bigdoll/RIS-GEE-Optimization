function G2a = G2func_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
% G2FUNC_ACTIVE2  Concave term F2(X).
%
%   function G2a = G2func_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (50)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    G2a = 0; % initialization
    N =  length(X);
    diag_X = X.*eye(N);
    W = sigma_sq*eye(NR) + sigma_sqris*G*diag_X*G';
    
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

        G2a = G2a + BW*log(real(det(W + pakm_sum)))/log(2); 
        
    end
   
end
