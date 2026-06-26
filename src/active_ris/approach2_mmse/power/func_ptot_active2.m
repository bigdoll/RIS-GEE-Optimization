function ptot_active2 = func_ptot_active2(p, X, H, sigma_sqris, K, Pca, mu)
% FUNC_PTOT_ACTIVE2  Total consumed power as a function of p (MMSE branch).
%
%   function ptot_active2 = func_ptot_active2(p, X, H, sigma_sqris, K, Pca, mu)
%
% Paper reference: Eqs. (10)-(11)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
        psum = 0;
        N = length(X);
        for k=1:K
            hk = H(:,k); % channel for user m b/w UE-RIS
            Hk = diag(hk); % diagonal matrix for hm
            mu_keq = mu + real(trace(Hk*(X - eye(N))*Hk')); 
            psum =  psum + mu_keq*p(k);
        end
        Pceq = sigma_sqris*real(trace(X - eye(N))) + Pca;
        ptot_active2 = psum + Pceq;
end
