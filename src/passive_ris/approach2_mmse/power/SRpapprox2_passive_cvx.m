function srp2_approx_cvx = SRpapprox2_passive_cvx(p, p_bar, G, H, X, sigma_sq, K, NR, BW, option)
% SRPAPPROX2_PASSIVE_CVX  CVX-expression form of the power surrogate (MMSE branch).
%
%   function srp2_approx_cvx = SRpapprox2_passive_cvx(p, p_bar, G, H, X, sigma_sq, K, NR, BW, option)
%
% Paper reference: Eq. (55)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
   
    f1 = func2_g1(p, G, H, X, sigma_sq, K, NR, BW, option);
    grad_f2 = func2_gradient(p_bar, G, H, X, sigma_sq, K, NR, BW);
    
    srp2_approx_cvx = f1 - grad_f2'* p ;
    
end
