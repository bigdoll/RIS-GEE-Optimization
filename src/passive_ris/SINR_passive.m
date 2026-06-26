function sinr_p = SINR_passive(p, C, G, H, gamma, sigma_sq, K)
% SINR_PASSIVE  Per-user SINR with a nearly-passive RIS (no RIS noise amplification).
%
%   function sinr_p = SINR_passive(p, C, G, H, gamma, sigma_sq, K)
%
% Paper reference: Eq. (2), sigma_RIS=0
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    %initialization
    sinr_p = zeros(K,1);
    
    for k=1:K
        hk = H(:,k); % channel for user k b/w UE-RIS
        Hk = diag(hk); % diagonal matrix for hk
        Ak = G * Hk; % Cascaded channel matrix for user k
        num_k = p(k) * abs(C(:,k)'* Ak * gamma)^2; % useful signal power for user k
        interf_m = 0; % initializing the interference + noise ratio for user k 
        for m = 1:K
            if m ~= k
                hm = H(:,m);
                Hm = diag(hm);
                Am = G * Hm;
                interf_m = interf_m + (p(m) * abs(C(:,k)'* Am * gamma)^2); 
                %continue % skip when m == k (exclude user k's signal power
            end     
        end
        denom_k =  sigma_sq*norm(C(:,k))^2 + interf_m; % interference + noise ratio for user k 
        sinr_p(k) = num_k / denom_k; % passive sinr for user k
    end
end
