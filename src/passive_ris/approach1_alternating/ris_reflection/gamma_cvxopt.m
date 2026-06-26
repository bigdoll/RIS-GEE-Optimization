function [gamma_sol, sr_sol, gee_sol, ng_iter] = gamma_cvxopt(p, G, H, gamma, sigma_sq, K, NR, PR, BW, mu, Pc, cons_mode)
% GAMMA_CVXOPT  ALGORITHM 1 (nearly-passive): sequential fractional-programming RIS optimization.
%
%   function [gamma_sol, sr_sol, gee_sol, ng_iter] = gamma_cvxopt(p, G, H, gamma, sigma_sq, K, NR, PR, BW, mu, Pc, cons_mode)
%
% cons_mode (optional): 'global' (default) enforces the global reflection
%   constraint  ||gamma||^2 <= N*PR  (a single budget over all elements);
%   'local' enforces the classical per-element constraint  |gamma_n| <= 1.
%   The 'local' mode is used to draw the comparison curve of Fig. 6.
%
% Paper reference: Algorithm 1, Problem (38)  (global vs. local: Sec. III, Fig. 6)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    if nargin < 12 || isempty(cons_mode); cons_mode = 'global'; end

    % Initialization
    ng_iter = 0;
    srapp_prev = 0;
    N = length(gamma);
    gamma_nxt = gamma;
    status = 'Solved';
    cvx_status = 'Failed';
    srapp_nxt = SRgapprox_passive(p, G, H, gamma_nxt, gamma_nxt, sigma_sq, K, NR, BW);

    %tolerance
    tol = 1e-6;
    
    while srapp_nxt - srapp_prev > tol || isempty(regexpi(cvx_status, 'Solved'))
    
        ng_iter = ng_iter + 1;
        srapp_prev = srapp_nxt;
        gamma_prev = gamma_nxt;
            
        cvx_begin
            cvx_quiet(true)
            cvx_precision high
            cvx_solver mosek
            variable gamma_opt(N) complex; % semidefinite % or hermitian semidefinite

            objective = SRgapprox_passive_cvx(p, G, H, gamma_opt, gamma_prev, sigma_sq, K, NR, BW);
            
            maximize objective;
            subject to
                if strcmpi(cons_mode, 'local')
                    abs(gamma_opt) <= 1;              % per-element |gamma_n| <= 1
                else
                    sum_square_abs(gamma_opt) <= N*PR; % global budget ||gamma||^2 <= N*PR
                end

        cvx_end
        
        % Robustness guard: keep the last good iterate on a failed/non-finite solve.
        if isempty(regexpi(cvx_status, 'Solved')) || any(~isfinite(gamma_opt(:)))
            warning('gamma_cvxopt:solveFallback', 'CVX returned "%s" at iteration %d; keeping previous iterate.', cvx_status, ng_iter);
            gamma_nxt = gamma_prev;
            cvx_status = 'Solved';
            break;
        end

        gamma_nxt = gamma_opt;
        srapp_nxt = SRgapprox_passive(p, G, H, gamma_nxt, gamma_prev, sigma_sq, K, NR, BW); %cvxopt_val 
    
    end   
    
    sr_sol = srapp_nxt;
    gamma_sol = gamma_nxt;
    
    if srapp_nxt < srapp_prev
        sr_sol = srapp_prev;
        gamma_sol = gamma_prev;
    end
    
    gee_sol = sr_sol / func_ptot(p, K, Pc, mu);
    
    fprintf('Gamma Opt: Number of Steps: %d, SR_opt: %f, GEE _opt: %f, Status of Solver: %s\n',ng_iter, sr_sol, gee_sol, cvx_status)
end
