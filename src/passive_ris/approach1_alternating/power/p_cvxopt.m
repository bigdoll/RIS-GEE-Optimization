function [p_sol, sr_sol, gee_sol, np_iter] = p_cvxopt(p, pmax, G, H, gamma, sigma_sq, K, NR, BW, mu, Pc, opt_bool)
% P_CVXOPT  ALGORITHM 2 (nearly-passive): sequential power optimization.
%
%   function [p_sol, sr_sol, gee_sol, np_iter] = p_cvxopt(p, pmax, G, H, gamma, sigma_sq, K, NR, BW, mu, Pc, opt_bool)
%
% Paper reference: Algorithm 2, Problem (43)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    mu0 = mu;
    Pc0 = Pc;
    
    if opt_bool == 0
        mu = 0;
        Pc = 1;
    end

    % Initialization
    np_iter = 0;
    p_sol = p;
    status = 'Solved';
    cvx_status = 'Failed';
    ln_prev = 0;
    sr_nxt = SRpapprox_passive(p_sol, p_sol, G, H, gamma, sigma_sq, K, NR, BW);
    ln_nxt = sr_nxt / func_ptot(p_sol, K, Pc, mu);

    %tolerance
    tol = 1e-6;
     
    while ln_nxt - ln_prev > tol || isempty(regexpi(cvx_status, 'Solved'))
    
        np_iter = np_iter + 1;
        sr_prev = sr_nxt;
        ln_prev = ln_nxt; 
        p_prev = p_sol;

        cvx_begin
            cvx_quiet(true)
            cvx_precision high
            cvx_solver mosek
            variable p_opt(K); % semidefinite % or hermitian semidefinite
            
            objective = SRpapprox_passive_cvx(p_opt, p_prev, G, H, gamma, sigma_sq, K, NR, BW) - opt_bool*ln_prev*(sum(mu*p_opt) + Pc);
           
            maximize objective;
            subject to
                sum(p_opt) <= pmax;
                p_opt >= 0;

        cvx_end
        
        % Robustness guard: keep the last good iterate on a failed/non-finite solve.
        if isempty(regexpi(cvx_status, 'Solved')) || any(~isfinite(p_opt(:)))
            warning('p_cvxopt:solveFallback', 'CVX returned "%s" at iteration %d; keeping previous iterate.', cvx_status, np_iter);
            p_sol = p_prev;
            cvx_status = 'Solved';
            break;
        end

        p_sol = abs(p_opt);
        sr_nxt = SRpapprox_passive(p_sol, p_prev, G, H, gamma, sigma_sq, K, NR, BW); %cvx_optval;
        ln_nxt = sr_nxt / func_ptot(p_sol, K, Pc, mu);
             
    end
    
    sr_sol = sr_nxt;
    gee_sol = ln_nxt;

    if ln_nxt < ln_prev
        fprintf('True, ln_nxt is less than ln_prev \n')
        sr_sol = sr_prev;
        gee_sol = ln_prev;
        p_sol = p_prev;
    end  
    
    if opt_bool == 0
        gee_sol = sr_sol / func_ptot(p_sol, K, Pc0, mu0);
    end
    
    fprintf('1st Approach -> Passive Power Opt: Number of Steps: %d, SR_opt: %f, GEE_opt: %f, Status of Solver: %s\n', np_iter, sr_sol, gee_sol, cvx_status)
end
