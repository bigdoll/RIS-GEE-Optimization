function grad_g2a = gradg2_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
% GRADG2_ACTIVE2  Gradient of f2(p) (MMSE branch) for the linearization.
%
%   function grad_g2a = gradg2_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (55)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    grad_g2a = zeros(K,1);
    gamma_vec = diag(X, 0);
    Gamma_herm = diag(gamma_vec); % Gamma_herm = Gamma*Gamma'
    W = sigma_sq*eye(NR) + sigma_sqris*G*Gamma_herm*G';

    for i=1:K
        hi = H(:,i); % channel for user i b/w UE-RIS
        Hi = diag(hi); % diagonal matrix for hi
        Ai = G * Hi; % Cascaded channel matrix for user i
        for k=1:K
            if k ~= i
                pakm_sum = 0;
                for m=1:K
                    if m ~= k
                        hm = H(:,m); % channel for user m b/w UE-RIS
                        Hm = diag(hm); % diagonal matrix for hm
                        Am = G * Hm; % Cascaded channel matrix for user m
                        akm = Am*X*Am';
                        pm = p(m);
                        pakm_sum = pakm_sum + pm*akm;
                    end
                end
                G2_inner = W + pakm_sum;
                grad_g2a(i) = grad_g2a(i) + real(trace(G2_inner\(Ai*X*Ai')));
            end
        end
        grad_g2a(i) = (BW/log(2)) * grad_g2a(i);
    end
end
