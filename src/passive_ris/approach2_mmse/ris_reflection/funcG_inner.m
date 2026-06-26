function [G1_inner, G2_inner] = funcG_inner(p, G, H, X, sigma_sq, K, NR, m)
% FUNCG_INNER  Inner per-user log-det term of the MMSE sum-rate.
%
%   function [G1_inner, G2_inner] = funcG_inner(p, G, H, X, sigma_sq, K, NR, m)
%
% Paper reference: Eq. (50)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    %initialization
    XG1_in = zeros(NR,NR,K);
    XG2_in = zeros(NR,NR,K);
    for k = 1:K
        hk = H(:,k); % channel for user k b/w UE-RIS
        Hk = diag(hk); % diagonal matrix for hk
        Ak = G * Hk; % Cascaded channel matrix
        pk = p(k);
        XG1_in(:,:,k) = pk*Ak*X*Ak';
        XG2_in(:,:,k) = XG1_in(:,:,k);
        if k == m
             XG2_in(:,:,k) = 0;
        end
    end
        G1_inner = sigma_sq*eye(NR) + sum(XG1_in,3);
        G2_inner = sigma_sq*eye(NR) + sum(XG2_in,3);
end
