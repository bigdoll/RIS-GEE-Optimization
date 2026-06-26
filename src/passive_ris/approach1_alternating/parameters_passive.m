function [A_bar, B_bar, D_bar, E_bar, F_bar] = parameters_passive(p, C, G, H, gamma_bar, sigma_sq, K)
% PARAMETERS_PASSIVE  Sequential-approximation constants A_bar..F_bar for the RIS subproblem.
%
%   function [A_bar, B_bar, D_bar, E_bar, F_bar] = parameters_passive(p, C, G, H, gamma_bar, sigma_sq, K)
%
% Paper reference: Eqs. (34)-(36)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    % parameters definition
    BW0 = 1;
    sinr_p = SINR_passive(p, C, G, H, gamma_bar, sigma_sq, K);
    [A_bar, ~] = data_rate(sinr_p, K, BW0);
    
    B_bar = BW0 * sinr_p;
    
    D_bar = zeros(K,1);
    for k=1:K
        ck = C(:,k);
        hk = H(:,k); % channel for user k b/w UE-RIS
        Hk = diag(hk); % diagonal matrix for hk
        Ak = G * Hk; % Cascaded channel matrix
        D_bar(k) = 2 / abs(ck'*Ak*gamma_bar);
    end
    
    E_bar = zeros(K,1);
    for k=1:K
        ck = C(:,k);
        pakm_sum = 0;
        for m = 1:K
            hm = H(:,m); % channel for user k b/w UE-RIS
            Hm = diag(hm); % diagonal matrix for hk
            Am = G * Hm; % Cascaded channel matrix
            akm = abs(ck'*Am*gamma_bar)^2;
            pakm_sum = pakm_sum + p(m)*akm;
        end
        E_bar(k) = 1 / (sigma_sq*(norm(ck))^2 + pakm_sum);
    end
    
    
    F_bar = zeros(K,1);
    for k=1:K
        ck = C(:,k);
        F_bar(k) = E_bar(k) * sigma_sq * (norm(ck))^2 + 1;
    end
    
end
