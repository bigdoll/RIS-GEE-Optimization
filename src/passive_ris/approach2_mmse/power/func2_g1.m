function g1 = func2_g1(p, G, H, X, sigma_sq, K, NR, BW, option)
% FUNC2_G1  Concave term g1(p) of the MMSE-embedded sum-rate.
%
%   function g1 = func2_g1(p, G, H, X, sigma_sq, K, NR, BW, option)
%
% Paper reference: Eq. (55)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    %initialization
    g1 = 0;
    
    for k = 1:K
        g1_inner = zeros(NR,NR);
        for m = 1:K
            hm = H(:,m); % channel for user k b/w UE-RIS
            Hm = diag(hm); % diagonal matrix for hk
            Am = G * Hm; % Cascaded channel matrix
            pm = p(m)/sigma_sq;
            g1_inner = g1_inner + pm*Am*X*Am';
        end
        
        g1_inner = 0.5*(g1_inner + g1_inner');
        
        if option == 1
            Xg1 = eye(NR) + g1_inner;
            g1 = g1 + BW*log_det(Xg1)/log(2); % BW*log(sigma_sq)/log(2)
        else
            g1_det = real(det(sigma_sq*(eye(NR) + g1_inner)));
            g1 = g1 + BW*log2(g1_det);
        end
    end
end
