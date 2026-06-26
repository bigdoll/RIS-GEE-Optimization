function Xtot_active = func_Xtot_active2(p, X, H, sigma_sqris, K, Pca, mu)
% FUNC_XTOT_ACTIVE2  Total consumed power as a function of X:  tr(R X) + ... .
%
%   function Xtot_active = func_Xtot_active2(p, X, H, sigma_sqris, K, Pca, mu)
%
% Paper reference: Eq. (51)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
        gsum = 0;
        N =  length(X);
        R = func_R(p, H, sigma_sqris, K, N);
        for k=1:K
            hk = H(:,k); % channel for user m b/w UE-RIS
            gsum =  gsum + p(k)*(mu - norm(hk)^2);
        end
        Xtot_active = real(trace(R*X)) + gsum + Pca - sigma_sqris*N;
end
