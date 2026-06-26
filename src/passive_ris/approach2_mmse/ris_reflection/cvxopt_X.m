function [X_sol, sr_sol, gee_sol, nX_iter] = cvxopt_X(p, G, H, X, sigma_sq, K, NR, PR, BW, mu, Pc)
% CVXOPT_X  ALGORITHM 4 (nearly-passive): semidefinite-relaxation RIS optimization in X.
%
%   function [X_sol, sr_sol, gee_sol, nX_iter] = cvxopt_X(p, G, H, X, sigma_sq, K, NR, PR, BW, mu, Pc)
%
% Paper reference: Algorithm 4, Problem (54)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    % Initialization
    nX_iter = 0;
    srapp_prev = 0;
    X_nxt = X;
    N =  length(X);
    status = 'Solved';
    cvx_status = 'Failed';

    srapp_nxt = func_SRmmse_approx(p, G, H, X_nxt, X_nxt, sigma_sq, K, NR, BW);

    %tolerance
    tol = 1e-6;
    
    while srapp_nxt - srapp_prev > tol || isempty(regexpi(cvx_status, 'Solved')) 
    
        nX_iter = nX_iter + 1;
        srapp_prev = srapp_nxt;
        X_prev = X_nxt;
        G2_grad = funcG2_grad(p, G, H, X_prev, sigma_sq, K, NR, BW);

        cvx_begin
            cvx_quiet(true)
            cvx_precision high
            cvx_solver mosek
            variable X_opt(N,N) hermitian semidefinite %complex semidefinite; % or hermitian semidefinite
            
            objective = funcG1_cvx(p, G, H, X_opt, sigma_sq, K, NR, BW) - real(trace(G2_grad'*X_opt));
      
            maximize objective;
            subject to
                trace(X_opt) <= N*PR;

        cvx_end
        
        % Robustness guard: keep the last good iterate on a failed/non-finite solve.
        if isempty(regexpi(cvx_status, 'Solved')) || any(~isfinite(X_opt(:)))
            warning('cvxopt_X:solveFallback', 'CVX returned "%s" at iteration %d; keeping previous iterate.', cvx_status, nX_iter);
            X_nxt = X_prev;
            cvx_status = 'Solved';
            break;
        end

        X_nxt = X_opt;
        srapp_nxt = func_SRmmse_approx(p, G, H, X_nxt, X_prev, sigma_sq, K, NR, BW);

    end
    
    X_sol = X_nxt;
    sr_sol = srapp_nxt;

    if srapp_nxt < srapp_prev
        X_sol = X_prev;
        sr_sol = srapp_prev;
    end
    
    gee_sol = sr_sol / func_ptot(p, K, Pc, mu);
    
    fprintf('2nd Approach -> Passive Xopt: Number of Steps: %d, SR_opt: %f, GEE_opt: %f, Status of Solver: %s\n', nX_iter, sr_sol, gee_sol, cvx_status)
end
