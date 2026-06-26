function [X_sol, srapp_sol, ln_sol, ng_iter] = X_cvxopt_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, PR, BW, mu, Pca, opt_bool)
% X_CVXOPT_ACTIVE2  ALGORITHM 4 (active): semidefinite-relaxation RIS optimization in X.
%
%   function [X_sol, srapp_sol, ln_sol, ng_iter] = X_cvxopt_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, PR, BW, mu, Pca, opt_bool)
%
% Paper reference: Algorithm 4, Problem (54)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    
    % Initialization
    ng_iter = 0;
    X_sol = X;
    N = length(X);
    status = 'Solved';
    cvx_status = 'Failed';
    srapp_prev = 0;
    ln_prev = 0;
    srapp_nxt = SRXapprox_active2(p, G, H, X, X, sigma_sq, sigma_sqris, K, NR, BW);
    ln_nxt = srapp_nxt / func_Xtot_active2(p, X, H, sigma_sqris, K, Pca, mu);
    X_prev =  X;
    
     if opt_bool == 0
        ln_nxt =  srapp_nxt;    
    end
    
    % computing Rnorm and PRmax
    R = func_R(p, H, sigma_sqris, K, N);
    rmax = abs(max(R, [], 'all'));
    Rnorm = R / rmax;
    n = 4;
    PRmax = sqrt((1/n)*trace(Rnorm) * N * PR);
    
    %tolerance
    tol = 1e-6;
    
    while ln_nxt - ln_prev > tol || isempty(regexpi(cvx_status, 'Solved'))
    
        ng_iter = ng_iter + 1;
        srapp_prev = srapp_nxt;
        ln_prev = ln_nxt;
        X_prev = X_sol;

        grad_G2 = G2grad_active2(p, G, H, X_prev, sigma_sq, sigma_sqris, K, NR, BW);
        
        cvx_begin
            cvx_quiet(true)
            cvx_precision high
            cvx_solver mosek
            variable X_opt(N,N) hermitian semidefinite; % semidefinite % or hermitian semidefinite
            
            G1 = G1func_active_cvx2(p, G, H, X_opt, sigma_sq, sigma_sqris, K, NR, BW);
            objective = G1 - real(trace(grad_G2'*X_opt)) - opt_bool*ln_prev*func_Xtot_active2(p, X_opt, H, sigma_sqris, K, Pca, mu);
%             objective = SRXapprox_active_cvx2(p, G, H, X_opt, X_prev, sigma_sq, sigma_sqris, K, NR, BW) - opt_bool*ln_prev*func_Xtot_active2(p, X_opt, H, sigma_sqris, K, Pca, mu);
            
            maximize objective;
            subject to
                real(trace(Rnorm*X_opt)) >= trace(Rnorm);
                real(trace(Rnorm*X_opt)) <= PRmax;
                
        cvx_end
        
        
        % Robustness guard: keep the last good iterate on a failed/non-finite solve.
        if isempty(regexpi(cvx_status, 'Solved')) || any(~isfinite(X_opt(:)))
            warning('X_cvxopt_active2:solveFallback', 'CVX returned "%s" at iteration %d; keeping previous iterate.', cvx_status, ng_iter);
            X_sol = X_prev;
            cvx_status = 'Solved';
            break;
        end

        X_sol = X_opt;
        srapp_nxt = SRXapprox_active2(p, G, H, X_sol, X_prev, sigma_sq, sigma_sqris, K, NR, BW);
        ln_nxt = srapp_nxt / func_Xtot_active2(p, X_sol, H, sigma_sqris, K, Pca, mu);
        
        if opt_bool == 0
            ln_nxt =  srapp_nxt;    
        end
%         X_norm = trace(X_sol);
%         X_cons = real(trace(Rnorm*X_sol));
%         fprintf('--- Step: %d, X_norm: %f, X_cons: %f, srapp_nxt: %f, ln_nxt: %f, Status of Solver: %s ---\n',ng_iter, X_norm, X_cons, srapp_nxt, ln_nxt, cvx_status)
    end   
    
    srapp_sol = srapp_nxt;
    ln_sol = ln_nxt;

    if ln_nxt < ln_prev
        fprintf('True, ln_nxt is less than ln_prev \n')
        X_sol = X_prev;
        srapp_sol = srapp_prev;
        ln_sol = ln_prev;
    end
    
     if opt_bool == 0
        ln_sol =  srapp_sol / func_Xtot_active2(p, X_sol, H, sigma_sqris, K, Pca, mu);    
     end
    
    fprintf('2nd approach -> Active Xopt: Number of Steps: %d, SR_opt: %f, GEE_opt: %f, Status of Solver: %s\n',ng_iter, srapp_sol, ln_sol, cvx_status) %
end
