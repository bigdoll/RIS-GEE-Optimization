function grad_G2 = funcG2_grad(p, G, H, X, sigma_sq, K, NR, BW)
% FUNCG2_GRAD  Gradient grad F2(X) for the linearization.
%
%   function grad_G2 = funcG2_grad(p, G, H, X, sigma_sq, K, NR, BW)
%
% Paper reference: Eq. (53)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    %Initialization
    N = length(X);
    grad_G2 = zeros(N,N);

    for k = 1:K
        grad_inner = zeros(N,N);
        for m = 1:K
            if m ~= k
                hm = H(:,m); % channel for user m b/w UE-RIS
                Hm = diag(hm); % diagonal matrix for hm
                Am = G * Hm; % Cascaded channel matrix at m
                [~,G2_inner] = funcG_inner(p, G, H, X, sigma_sq, K, NR, k);
                C_inv = G2_inner\Am;
                pm = p(m); 
                grad_inner = grad_inner + pm*Am'*C_inv;
            end
        end
        grad_G2 = grad_G2 + BW*grad_inner/log(2);
    end
end
