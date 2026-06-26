function [srapprox_pcvx] = SRgapprox_passive_cvx(p, G, H, gamma, gamma_bar, sigma_sq, K, NR, BW)
% SRGAPPROX_PASSIVE_CVX  CVX-expression form of the gamma surrogate.
%
%   function [srapprox_pcvx] = SRgapprox_passive_cvx(p, G, H, gamma, gamma_bar, sigma_sq, K, NR, BW)
%
% Paper reference: Eqs. (35)-(36)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    %initialization
    srapprox_pcvx = 0;
    
    % Linear MMSE update
    C = LMMSE_receiver_passive(p, G, H, gamma_bar, sigma_sq, K, NR);
    
    %paramters
    [~, B_bar, D_bar, E_bar, ~] = parameters_passive(p, C, G, H, gamma_bar, sigma_sq, K);
    
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
        Zk = abs(ck'*Ak*gamma_bar);
        deriv = (Ak'*(ck*ck')*Ak)*gamma_bar / Zk;
        srapprox_pcvx = srapprox_pcvx + BW*B_bar(k)*(D_bar(k)*real(deriv'*gamma) - E_bar(k)*Yp);
    end
end
