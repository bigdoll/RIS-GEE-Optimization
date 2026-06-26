function [A_bar, B_bar, D_bar, E_bar, F_bar] = parameters_active(p, C, G, H, gamma_bar, sigma_sq, sigma_sqris, K)
% PARAMETERS_ACTIVE  Sequential-approximation constants A_bar,B_bar,D_bar,E_bar,F_bar for the RIS subproblem.
%
%   function [A_bar, B_bar, D_bar, E_bar, F_bar] = parameters_active(p, C, G, H, gamma_bar, sigma_sq, sigma_sqris, K)
%
% Paper reference: Eqs. (34)-(36)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    % parameters definition
    BW0 = 1;
    sinr_a = SINR_active(p, C, G, H, gamma_bar, sigma_sq, sigma_sqris, K);
    [A_bar, ~] = data_rate_active(sinr_a, K, BW0);
    
    B_bar = BW0 * sinr_a;
    
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
        uk = G'*ck;
        Uk_tilt = diag(abs(uk).^2);
        noise_k = sigma_sq*(norm(ck))^2 + sigma_sqris*real(gamma_bar'*Uk_tilt*gamma_bar);
        E_bar(k) = 1 / (noise_k + pakm_sum);
    end
    
    
    F_bar = zeros(K,1);
    for k=1:K
        ck = C(:,k);
        F_bar(k) = E_bar(k) * sigma_sq * (norm(ck))^2 + 1;
    end
    
end
