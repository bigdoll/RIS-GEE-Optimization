function [G1, G2] = func_G(p, G, H, X, sigma_sq, K, NR, BW)
% FUNC_G  Concave term F1(X) of the MMSE sum-rate.
%
%   function [G1, G2] = func_G(p, G, H, X, sigma_sq, K, NR, BW)
%
% Paper reference: Eq. (50)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    %Initialization
    G1 = 0;
    G2 = 0;
    
    for k = 1:K
       [G1_inner, G2_inner] = funcG_inner(p, G, H, X, sigma_sq, K, NR, k);   
       G1 = G1 + BW*log2(real(det(G1_inner)));
       G2 = G2 + BW*log2(real(det(G2_inner)));
    end
end
