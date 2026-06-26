function C_active = LMMSE_receiver_active(p, G, H, gamma, sigma_sq, sigma_sqris, K, NR)
% LMMSE_RECEIVER_ACTIVE  Closed-form linear MMSE receive filters c_k (active RIS).
%
%   function C_active = LMMSE_receiver_active(p, G, H, gamma, sigma_sq, sigma_sqris, K, NR)
%
% Paper reference: Eq. (44)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    % initializing the receive beamformer vectors
    C_active = zeros(NR,K);    
    Gamma = diag(gamma);
    W = sigma_sq*eye(NR) + sigma_sqris*G*(Gamma*Gamma')*G';

    for k=1:K
        hk = H(:,k); % channel for user k b/w UE-RIS
        Hk = diag(hk); % diagonal matrix for hk
        Ak = G * Hk; % Cascaded channel matrix
        C_inner = 0; %zeros(NR,NR); 
        for m=1:K
            if m ~= k
                hm = H(:,m);
                Hm = diag(hm);
                Am = G * Hm;
                C_inner = C_inner + (p(m) * Am * (gamma * gamma') * Am');
            end
        end
        
        Ckmat = C_inner + W;
        Ag = Ak * gamma;
        C_active(:,k) = sqrt(p(k)) * pinv(Ckmat)*Ag;
    end
    
end
