function [srapprox_acvx] = SRgapprox_active_cvx(p, G, H, gamma, gamma_bar, sigma_sq, sigma_sqris, K, NR, BW)
% SRGAPPROX_ACTIVE_CVX  CVX-expression form of the gamma surrogate (objective inside cvx_begin).
%
%   function [srapprox_acvx] = SRgapprox_active_cvx(p, G, H, gamma, gamma_bar, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eqs. (35)-(36)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    %initialization
    srapprox_acvx = 0;
    
    C = LMMSE_receiver_active(p, G, H, gamma_bar, sigma_sq, sigma_sqris, K, NR);
    
    %paramters
    [~, B_bar, D_bar, E_bar, ~] = parameters_active(p, C, G, H, gamma_bar, sigma_sq, sigma_sqris, K);
    
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
             pm = p(m);
             Yp = Yp + pm * square_abs(ck'*Am*gamma);
        end
      
        uk = G'*ck;
        Uk_tilt = diag(abs(uk).^2);
        Zk = abs(ck'*Ak*gamma_bar);
        deriv = (Ak'*(ck*ck')*Ak)*gamma_bar / Zk;
        srapprox_acvx = srapprox_acvx + BW*B_bar(k)*(D_bar(k)*real(deriv'*gamma) - E_bar(k)*(sigma_sqris*real(gamma'*Uk_tilt*gamma) + Yp));
    end
end
