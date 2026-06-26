function [srapprox_a] = SRgapprox_active(p, G, H, gamma, gamma_bar, sigma_sq, sigma_sqris, K, NR, BW)
% SRGAPPROX_ACTIVE  Concave surrogate of the MMSE sum-rate in gamma (true value).
%
%   function [srapprox_a] = SRgapprox_active(p, G, H, gamma, gamma_bar, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eqs. (35)-(36)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    %initialization
    srapprox_a = 0;

    % Linear MMSE update
    C = LMMSE_receiver_active(p, G, H, gamma_bar, sigma_sq, sigma_sqris, K, NR);
    
    %paramters
    [A_bar, B_bar, D_bar, E_bar, F_bar] = parameters_active(p, C, G, H, gamma_bar, sigma_sq, sigma_sqris, K);
    
    for k=1:K
        ck = C(:,k);
        hk = H(:,k); % channel for user k b/w UE-RIS
        Hk = diag(hk); % diagonal matrix for hk
        Ak = G * Hk; % Cascaded channel matrix
        Yp  = 0;
        for m = 1:K
             hm = H(:,m); % channel for user k b/w UE-RIS
             Hm = diag(hm); % diagonal matrix for hk
             Am = G * Hm; % Cascaded channel matrix
             Yp = Yp + p(m) * square_abs(ck'*Am*gamma);
        end
        uk = G'*ck;
        Uk_tilt = diag(abs(uk).^2);
        Zk = abs(ck'*Ak*gamma_bar);
        deriv = (Ak'*(ck*ck')*Ak)*gamma_bar / Zk;
        srapprox_a = srapprox_a + BW*(A_bar(k) + B_bar(k)*(D_bar(k)*(Zk + real(deriv'*(gamma-gamma_bar))) - E_bar(k)*(sigma_sqris*real(gamma'*Uk_tilt*gamma) + Yp) - F_bar(k)));
    end
end
