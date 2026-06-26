function gtot = func_gtot_active_cvx(p, gamma, H, sigma_sqris, K, Pca, mu)
% FUNC_GTOT_ACTIVE_CVX  CVX-expression form of the GEE denominator in gamma.
%
%   gtot = func_gtot_active_cvx(p, gamma, H, sigma_sqris, K, Pca, mu)
%
% Identical to FUNC_GTOT_ACTIVE, but the quadratic term gamma'*R*gamma is
% written with the CVX atom quad_form(gamma, R) so that the disciplined-convex
% parser certifies it as a convex function of the CVX variable gamma. Use this
% version ONLY inside a cvx_begin ... cvx_end block (i.e. when gamma is a CVX
% variable); use FUNC_GTOT_ACTIVE for numeric evaluation.
%
% Paper reference: Eq. (11)  (denominator of the GEE, gamma subproblem)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    N = length(gamma);
    R = func_R(p, H, sigma_sqris, K, N);   % positive-definite, Eq. (6)

    gsum = 0;
    for k = 1:K
        gsum = gsum + p(k)*(mu - norm(H(:,k))^2);
    end

    % quad_form(gamma, R) == real(gamma'*R*gamma), but DCP-certified convex.
    gtot = quad_form(gamma, R) + gsum + Pca - sigma_sqris*N;
end
