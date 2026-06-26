function grad_g2 = gradient_g2(p, C, G, H, gamma, sigma_sq, K, BW)
% GRADIENT_G2  Gradient of f2(p) (numerator of g2) for the linearization.
%
%   function grad_g2 = gradient_g2(p, C, G, H, gamma, sigma_sq, K, BW)
%
% Paper reference: Eq. (42)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    grad_g2 = zeros(K,1);
    
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
                        pakm_sum = pakm_sum + p(m)*akm;
                    end
                end
                dk = sigma_sq*(norm(ck)^2);
                grad_g2(i) = grad_g2(i) + aki / (dk + pakm_sum);
            end
        end
        grad_g2(i) = (BW/log(2)) * grad_g2(i);
    end
    
end
