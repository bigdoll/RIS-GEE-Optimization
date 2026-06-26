function grad_g2 = func2_gradient(p, G, H, X, sigma_sq, K, NR, BW)
% FUNC2_GRADIENT  Gradient of f2(p) (MMSE branch) for the linearization.
%
%   function grad_g2 = func2_gradient(p, G, H, X, sigma_sq, K, NR, BW)
%
% Paper reference: Eq. (55)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    %Initialization
    grad_g2 = zeros(K,1);

    for i = 1:K
        hi = H(:,i); % channel for user k b/w UE-RIS
        Hi = diag(hi); % diagonal matrix for hk
        Ai = G * Hi; % Cascaded channel matrix at i
        for k = 1:K
            if k ~= i
                C_in = zeros(NR,NR);
                for m = 1:K
                    if m ~= k
                        hm = H(:,m); % channel for user k b/w UE-RIS
                        Hm = diag(hm); % diagonal matrix for hk
                        Am = G * Hm; % Cascaded channel matrix
                        C_in = C_in + p(m)*Am*X*Am';
                    end
                end
                C_out = sigma_sq*eye(NR) + C_in;
                Xtr = real(trace(C_out\(Ai*X*Ai')));
                grad_g2(i) = grad_g2(i) + BW*Xtr/log(2);
            end
        end
    end
end
