function sinr_a = SINR_active(p, C, G, H, gamma, sigma_sq, sigma_sqris, K)
% SINR_ACTIVE  Per-user SINR with an active RIS (includes RIS noise amplification).
%
%   function sinr_a = SINR_active(p, C, G, H, gamma, sigma_sq, sigma_sqris, K)
%
% Paper reference: Eq. (2)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    %initialization
    sinr_a = zeros(K,1);
    
    for k=1:K
        ck = C(:,k);
        hk = H(:,k); % channel for user k b/w UE-RIS
        Hk = diag(hk); % diagonal matrix for hk
        Ak = G * Hk; % Cascaded channel matrix for user k
        num_k = p(k) * abs(ck'* Ak * gamma)^2; % useful signal power for user k
        interf_m = 0; % initializing the interference term for user k (that is all other users except k)
        for m = 1:K
            if m ~= k
                hm = H(:,m);
                Hm = diag(hm);
                Am = G * Hm;
                interf_m = interf_m + (p(m) * abs(ck'* Am * gamma)^2); 
            end     
        end
        uk = G'*ck;
        Uk_tilt = diag(abs(uk).^2);
        noise_k = sigma_sq*(norm(ck))^2 + sigma_sqris*real(gamma'*Uk_tilt*gamma); % noise term for user k recieved signal
        denom_k =  noise_k + interf_m; % interference + noise ratio for user k 
        sinr_a(k) = num_k / denom_k; % active sinr for user k
    end
end
