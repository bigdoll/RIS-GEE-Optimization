function G1_cvx = funcG1_cvx(p, G, H, X, sigma_sq, K, NR, BW)
% FUNCG1_CVX  CVX-expression form of F1(X).
%
%   function G1_cvx = funcG1_cvx(p, G, H, X, sigma_sq, K, NR, BW)
%
% Paper reference: Eq. (50)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    % Initialization
    G1_cvx = 0;
    
    for i = 1:K
        XG1_tot = zeros(NR,NR);
        for k = 1:K
            hk = H(:,k); % channel for user k b/w UE-RIS
            Hk = diag(hk); % diagonal matrix for hk
            Ak = G * Hk; % Cascaded channel matrix
            pk = p(k)/sigma_sq; 
            XG1_tot = XG1_tot + pk*Ak*X*Ak';
        end
        XG1_tot = 0.5*(XG1_tot + XG1_tot');
        G1_inner = eye(NR) + XG1_tot;
        G1_cvx  = G1_cvx + BW*log_det(G1_inner)/log(2);
    end
end
