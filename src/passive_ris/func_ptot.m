function ptot = func_ptot(p, K, Pc, mu)
% FUNC_PTOT  Total consumed power  sum_k mu*p_k + Pc (nearly-passive RIS).
%
%   function ptot = func_ptot(p, K, Pc, mu)
%
% Paper reference: Eqs. (10),(14)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
        psum = 0;
        for k=1:K
            psum =  psum + mu*p(k);
        end
        ptot = psum + Pc;
end
