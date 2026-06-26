function C_passive = LMMSE_receiver_passive(p, G, H, gamma, sigma_sq, K, NR)
% LMMSE_RECEIVER_PASSIVE  Closed-form linear MMSE receive filters c_k (nearly-passive RIS).
%
%   function C_passive = LMMSE_receiver_passive(p, G, H, gamma, sigma_sq, K, NR)
%
% Paper reference: Eq. (44)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    % initializing the receive beamformer vectors
    C_passive = zeros(NR,K); 

    for k=1:K
        hk = H(:,k); % channel for user k b/w UE-RIS
        Hk = diag(hk); % diagonal matrix for hk
        Ak = G * Hk; % Cascaded channel matrix
        C_inner = zeros(NR,NR); 
        for m=1:K
            if m ~= k
                hm = H(:,m);
                Hm = diag(hm);
                Am = G * Hm;
                C_inner = C_inner + (p(m) * Am * (gamma * gamma') * Am');
            end
        end
        Ck_mat = C_inner + sigma_sq * eye(NR);
        C_passive(:,k) = sqrt(p(k)) * pinv(Ck_mat)*(Ak * gamma); % LMMSE recieve filter vector for user k 
    end
    
end
