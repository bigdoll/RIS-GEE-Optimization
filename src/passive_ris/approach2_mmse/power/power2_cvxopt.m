function [p_sol, sr_sol, gee_sol, np_iter] = power2_cvxopt(p, pmax, G, H, X, sigma_sq, K, NR, BW, mu, Pc, opt_bool)
% POWER2_CVXOPT  ALGORITHM 5 (nearly-passive): sequential power optimization (MMSE branch).
%
%   function [p_sol, sr_sol, gee_sol, np_iter] = power2_cvxopt(p, pmax, G, H, X, sigma_sq, K, NR, BW, mu, Pc, opt_bool)
%
% Paper reference: Algorithm 5, Problem (56)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    mu0 = mu;
    Pc0 = Pc;
    
    if opt_bool == 0
        mu = 0;
        Pc = 1;
    end

    % Initialization
    np_iter = 0;
    srapp_prev = 0;
    p_nxt = p;
    status = 'Solved';
    cvx_status = 'Failed';
    option = [0,1];

    srapp_nxt = SRpapprox2_passive(p_nxt, p_nxt, G, H, X, sigma_sq, K, NR, BW, option(1));
    ln_prev = 0;
    ln_nxt = srapp_nxt / func_ptot(p_nxt, K, Pc, mu);


    %tolerance
    tol = 1e-6;
    
    while ln_nxt - ln_prev > tol || isempty(regexpi(cvx_status, 'Solved')) 
    
        np_iter = np_iter + 1;
        srapp_prev = srapp_nxt;
        ln_prev = ln_nxt;
        p_prev = p_nxt;
            
        cvx_begin
            cvx_quiet(true)
            cvx_precision high
            cvx_solver mosek
            variable p_opt(K); % semidefinite % or hermitian semidefinite

            objective = SRpapprox2_passive_cvx(p_opt, p_prev, G, H, X, sigma_sq, K, NR, BW, option(2)) - opt_bool*ln_prev*(mu*sum(p_opt) + Pc);

            maximize objective;
            subject to
                sum(p_opt) <= pmax;
                p_opt >= 0;

        cvx_end
        
        % Robustness guard: keep the last good iterate on a failed/non-finite solve.
        if isempty(regexpi(cvx_status, 'Solved')) || any(~isfinite(p_opt(:)))
            warning('power2_cvxopt:solveFallback', 'CVX returned "%s" at iteration %d; keeping previous iterate.', cvx_status, np_iter);
            p_nxt = p_prev;
            cvx_status = 'Solved';
            break;
        end

        p_nxt = abs(p_opt);
        srapp_nxt = SRpapprox2_passive(p_nxt, p_prev, G, H, X, sigma_sq, K, NR, BW, option(1)); %cvx_optval;
        ln_nxt = srapp_nxt / func_ptot(p_nxt, K, Pc, mu);
        
    end
    
    p_sol = p_nxt;
    sr_sol = srapp_nxt;
    gee_sol = ln_nxt;

    if ln_nxt < ln_prev
        p_sol = p_prev;
        sr_sol = srapp_prev;
        gee_sol = ln_prev;
    end
    
    if opt_bool == 0
        gee_sol = sr_sol / func_ptot(p_sol, K, Pc0, mu0);
    end
    
    fprintf('2nd Approach -> Passive Power Opt: Number of Steps: %d, SR_opt: %f, GEE_opt: %f, Status of Solver: %s\n',np_iter, sr_sol, gee_sol, cvx_status)
end
