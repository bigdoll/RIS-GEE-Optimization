function gtot_active = func_gtot_active(p, gamma, H, sigma_sqris, K, Pca, mu)
% FUNC_GTOT_ACTIVE  Total consumed power as a function of gamma (RIS subproblem).
%
%   function gtot_active = func_gtot_active(p, gamma, H, sigma_sqris, K, Pca, mu)
%
% Paper reference: Eq. (11)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
        gsum = 0;
        N =  length(gamma);
        R = func_R(p, H, sigma_sqris, K, N);
        for k=1:K
            hk = H(:,k); % channel for user m b/w UE-RIS
%             Hk = diag(hk);% diagonal matrix for hk
            gsum =  gsum + p(k)*(mu - norm(hk)^2);
%             R = R + (p(k)*(Hk'*Hk) + sigma_sqris*eye(N));
        end
%         test1 = real(gamma'*R*gamma);
%         test2 = gsum + Pca - sigma_sqris*N;
        gtot_active = real(gamma'*R*gamma) + gsum + Pca - sigma_sqris*N;
end
