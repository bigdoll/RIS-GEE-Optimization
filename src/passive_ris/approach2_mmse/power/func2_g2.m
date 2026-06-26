function g2 = func2_g2(p, G, H, X, sigma_sq, K, NR, BW, option)
% FUNC2_G2  Convex term g2(p) of the MMSE-embedded sum-rate.
%
%   function g2 = func2_g2(p, G, H, X, sigma_sq, K, NR, BW, option)
%
% Paper reference: Eq. (55)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    %initialization
    g2 = 0;

    for k = 1:K
        g2_inner = zeros(NR,NR);
        for m = 1:K
            if m ~= k
                hm = H(:,m); % channel for user k b/w UE-RIS
                Hm = diag(hm); % diagonal matrix for hk
                Am = G * Hm; % Cascaded channel matrix
                g2_inner = g2_inner + p(m)*Am*X*Am';
            end
        end
        
        g2_inner = 0.5*(g2_inner +g2_inner');
        
        if option == 1
            Xg2 = sigma_sq*eye(NR) + g2_inner;
            g2 = g2 + BW*log_det(Xg2)/log(2);
        else
            g2_det = real(det(sigma_sq*eye(NR) + g2_inner));
            g2 = g2 + BW*log2(g2_det);
        end
    end
end
