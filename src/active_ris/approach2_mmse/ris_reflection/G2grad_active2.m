function grad_G2a = G2grad_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
% G2GRAD_ACTIVE2  Gradient grad F2(X) for the linearization.
%
%   function grad_G2a = G2grad_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW)
%
% Paper reference: Eq. (53)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    % Initialization
    N =  length(X);
    grad_G2a = zeros(N,N);
    diag_X = X.*eye(N);
    W = sigma_sq*eye(NR) + sigma_sqris*(G*diag_X*G');

    for k=1:K
        
        paki_inner = 0;
        for i=1:K
            if i ~= k
                hi = H(:,i); % channel for user m b/w UE-RIS
                Hi = diag(hi); % diagonal matrix for hm
                Ai = G * Hi; % Cascaded channel matrix for user m
                AXi = Ai*X*Ai';
                pi = p(i);
                paki_inner = paki_inner + pi*AXi;
            end
        end
        
        Tk = W + paki_inner;
       
       pakm_outer = 0;
       for m=1:K  
            if m ~= k
                hm = H(:,m); % channel for user m b/w UE-RIS
                Hm = diag(hm); % diagonal matrix for hm
                Am = G * Hm; % Cascaded channel matrix for user m
                ATm = Am'*(Tk\Am);
                pm = p(m);
                pakm_outer = pakm_outer + pm*ATm;
            end
       end
        
%         test = (G'*(Tk\G));
        GT_inner = sigma_sqris*(G'*(Tk\G)).*eye(N);
        grad_G2a = grad_G2a + (BW/log(2))*(GT_inner + pakm_outer);

    end
end
