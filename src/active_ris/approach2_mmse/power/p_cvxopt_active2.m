function [p_sol, srapp_sol, ln_sol, np_iter] = p_cvxopt_active2(p, pmax, G, H, X, sigma_sq, sigma_sqris, K, NR, BW, mu, Pca, opt_bool)
% P_CVXOPT_ACTIVE2  ALGORITHM 5 (active): sequential power optimization (MMSE branch).
%
%   function [p_sol, srapp_sol, ln_sol, np_iter] = p_cvxopt_active2(p, pmax, G, H, X, sigma_sq, sigma_sqris, K, NR, BW, mu, Pca, opt_bool)
%
% Paper reference: Algorithm 5, Problem (56)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    % Initialization
    np_iter = 0;
    srapp_prev = 0;
    p_sol = p;
    ln_prev = 0;
    srapp_nxt = SRpapprox_active2(p, p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW);
    ln_nxt = srapp_nxt / func_ptot_active2(p, X, H, sigma_sqris, K, Pca, mu);
    cvx_status = 'Failed';
    status = 'Solved';
    
     if opt_bool == 0
        ln_nxt =  srapp_nxt;    
    end

    %tolerance
    tol = 1e-6;
     
    while ln_nxt - ln_prev > tol || isempty(regexpi(cvx_status, 'Solved'))
    
        np_iter = np_iter + 1;
        srapp_prev = srapp_nxt;
        ln_prev = ln_nxt; 
        p_prev = p_sol;

        cvx_begin
            cvx_quiet(true)
            cvx_precision high
            cvx_solver mosek
            variable p_opt(K); % semidefinite % or hermitian semidefinite
            
            objective = SRpapprox_active_cvx2(p_opt, p_prev, G, H, X, sigma_sq, sigma_sqris, K, NR, BW) - opt_bool*ln_prev*func_ptot_active2(p_opt, X, H, sigma_sqris, K, Pca, mu);
            
            maximize objective;
            subject to
                sum(p_opt) <= pmax;
                p_opt >= 0;

        cvx_end
        
        % Robustness guard: keep the last good iterate on a failed/non-finite solve.
        if isempty(regexpi(cvx_status, 'Solved')) || any(~isfinite(p_opt(:)))
            warning('p_cvxopt_active2:solveFallback', 'CVX returned "%s" at iteration %d; keeping previous iterate.', cvx_status, np_iter);
            p_sol = p_prev;
            cvx_status = 'Solved';
            break;
        end

        p_sol = abs(p_opt);
        srapp_nxt = SRpapprox_active2(p_sol, p_prev, G, H, X, sigma_sq, sigma_sqris, K, NR, BW);
        ln_nxt = srapp_nxt / func_ptot_active2(p_sol, X, H, sigma_sqris, K, Pca, mu);
        
         if opt_bool == 0
            ln_nxt =  srapp_nxt;    
         end
             
    end
    
    srapp_sol = srapp_nxt;
    ln_sol = ln_nxt;

    if ln_nxt < ln_prev
        fprintf('True, ln_nxt is less than ln_prev \n')
        srapp_sol = srapp_prev;
        ln_sol = ln_prev;
        p_sol = p_prev;

    end  
    
    if opt_bool == 0
        ln_sol = srapp_sol / func_ptot_active2(p_sol, X, H, sigma_sqris, K, Pca, mu);
    end
     
    fprintf('2nd Approach -> Power Opt: Number of Steps: %d, SR_opt: %f, GEE_opt Value : %f, Status of Solver: %s\n',np_iter, srapp_sol, ln_sol, cvx_status)
end
