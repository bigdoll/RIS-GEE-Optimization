function ptot_active = func_ptot_active(p, gamma, H, sigma_sqris, K, Pca, mu)
% FUNC_PTOT_ACTIVE  Total consumed power as a function of p (power subproblem).
%
%   function ptot_active = func_ptot_active(p, gamma, H, sigma_sqris, K, Pca, mu)
%
% Paper reference: Eqs. (10)-(11),(39)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
        psum = 0;
        N = length(gamma);
        for k=1:K
            hk = H(:,k); % channel for user m b/w UE-RIS
            Hk = diag(hk); % diagonal matrix for hm
            mu_keq = mu + real(trace(Hk*((gamma*gamma') - eye(N))*Hk')); 
            psum =  psum + mu_keq*p(k);
        end
        Pceq = sigma_sqris*real(trace((gamma*gamma') - eye(N))) + Pca;
        ptot_active = psum + Pceq;
end
