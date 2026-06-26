function [p_sol, srapp_sol, ln_sol, np_iter] = p_cvxopt_active(p, pmax, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW, mu, Pca, opt_bool)
% P_CVXOPT_ACTIVE  ALGORITHM 2 (active): sequential fractional-programming power optimization.
%
%   function [p_sol, srapp_sol, ln_sol, np_iter] = p_cvxopt_active(p, pmax, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW, mu, Pca, opt_bool)
%
% Paper reference: Algorithm 2, Problem (43)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
%     
%     if opt_bool == 0
%         mu = 0;
%         Pca = 1;
%     end

    % Initialization
    np_iter = 0;
    srapp_prev = 0;
    p_sol = p;
    status = 'Solved';
    cvx_status = 'Failed';
    ln_prev = 0;
    srapp_nxt = SRpapprox_active(p, p, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW);
    ln_nxt = srapp_nxt / func_ptot_active(p, gamma, H, sigma_sqris, K, Pca, mu);
    
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
            
            objective = SRpapprox_active_cvx(p_opt, p_prev, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW) - opt_bool*ln_prev*func_ptot_active(p_opt, gamma, H, sigma_sqris, K, Pca, mu);
            % objective_old = g1_cvx(p_opt, C, G, H, gamma, sigma_sq, K, BW) - M'*p_opt;

            maximize objective;
            subject to
                sum(p_opt) <= pmax;
                p_opt >= 0;

        cvx_end

        % Robustness guard: keep the last good iterate on a failed/non-finite solve.
        if isempty(regexpi(cvx_status, 'Solved')) || any(~isfinite(p_opt))
            warning('p_cvxopt_active:solveFallback', 'CVX returned "%s" at iteration %d; keeping previous iterate.', cvx_status, np_iter);
            p_sol = p_prev;
            cvx_status = 'Solved';
            break;
        end

        p_sol = abs(p_opt);
        srapp_nxt = SRpapprox_active(p_sol, p_prev, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW); %cvx_optval;
        ln_nxt = srapp_nxt / func_ptot_active(p_sol, gamma, H, sigma_sqris, K, Pca, mu);
        
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
        ln_sol = srapp_sol / func_ptot_active(p_sol, gamma, H, sigma_sqris, K, Pca, mu);
    end

    fprintf('1st Approach -> Power Opt: Number of Steps: %d, GEE_opt Value : %d, Status of Solver: %s\n',np_iter, ln_sol, cvx_status)
end
