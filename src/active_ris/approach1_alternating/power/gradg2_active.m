function grad_g2a = gradg2_active(p, C, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW)
% GRADG2_ACTIVE  Gradient of f2(p) (numerator of g2) used for the linearization.
%
%   function grad_g2a = gradg2_active(p, C, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (42)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    grad_g2a = zeros(K,1);
    Gamma = diag(gamma);
    W = sigma_sq*eye(NR) + sigma_sqris*G*(Gamma*Gamma')*G';

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
                        pm = p(m);
                        pakm_sum = pakm_sum + pm*akm;
                    end
                end
                dk = real(ck'*W*ck);
                grad_g2a(i) = grad_g2a(i) + aki / (dk + pakm_sum);
            end
        end
        grad_g2a(i) = (BW/log(2)) * grad_g2a(i);
    end

end
