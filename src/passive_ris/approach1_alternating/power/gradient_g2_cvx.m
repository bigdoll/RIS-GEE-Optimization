function grad_g2_cvx = gradient_g2_cvx(p, C, G, H, gamma, sigma_sq, K, BW)
% GRADIENT_G2_CVX  CVX-expression form of grad f2(p).
%
%   function grad_g2_cvx = gradient_g2_cvx(p, C, G, H, gamma, sigma_sq, K, BW)
%
% Paper reference: Eq. (42)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    grad_g2_cvx = zeros(K,1);
    
    for i=1:K
        hi = H(:,i); % channel for user i b/w UE-RIS
        Hi = diag(hi); % diagonal matrix for hi
        Ai = G * Hi; % Cascaded channel matrix for user i
        for k=1:K
            if k ~= i
                ck = C(:,k); % Linear MMSE received filter for user k 
                aki = abs(ck'*Ai*gamma)^2;
                pakm_sum = 0;
                for m=1:K
                    if m ~= k
                        hm = H(:,m); % channel for user m b/w UE-RIS
                        Hm = diag(hm); % diagonal matrix for hm
                        Am = G * Hm; % Cascaded channel matrix for user m
                        akm = abs(ck'*Am*gamma)^2;
                        pm = p(m) / sigma_sq;
                        pakm_sum = pakm_sum + pm*akm;
                    end
                end
                dk = (norm(ck))^2;
                grad_g2_cvx(i) = grad_g2_cvx(i) + aki / (dk + pakm_sum);
            end
        end
        grad_g2_cvx(i) = (BW/(sigma_sq*log(2))) * grad_g2_cvx(i);
    end

end
