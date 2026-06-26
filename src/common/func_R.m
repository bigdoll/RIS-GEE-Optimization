function R =  func_R(p, H, sigma_sqris, K, N)
% FUNC_R  Positive-definite matrix  R = sum_k p_k H_k^H H_k + sigma_RIS^2 I_N.
%
%   function R =  func_R(p, H, sigma_sqris, K, N)
%
% Paper reference: Eq. (6)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    R = sigma_sqris * eye(N);
    for k =1:K
        hk = H(:,k);
        Hk = diag(hk);
        R = R + p(k)*(Hk'*Hk);
    end
end
