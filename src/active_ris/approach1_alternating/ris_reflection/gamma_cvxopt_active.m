function [gamma_sol, srapp_sol, ln_sol, ng_iter] = gamma_cvxopt_active(p, G, H, gamma, sigma_sq, sigma_sqris, K, NR, PR, BW, mu, Pca, opt_bool, cons_mode)
% GAMMA_CVXOPT_ACTIVE  ALGORITHM 1 (active): sequential fractional-programming RIS optimization.
%
%   function [gamma_sol, srapp_sol, ln_sol, ng_iter] = gamma_cvxopt_active(p, G, H, gamma, sigma_sq, sigma_sqris, K, NR, PR, BW, mu, Pca, opt_bool, cons_mode)
%
% cons_mode (optional): 'global' (default) enforces the global reflection
%   constraint  tr(R) <= tr(R*gamma*gamma') <= PRmax + tr(R)  (Eqs. 7-8);
%   'local' enforces the classical per-element constraint  |gamma_n| <= 1,
%   used to draw the comparison curve of Fig. 6.
%
% Paper reference: Algorithm 1, Problem (38)  (global vs. local: Sec. III, Fig. 6)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    if nargin < 14 || isempty(cons_mode); cons_mode = 'global'; end

    % Initialization
    ng_iter = 0;
    N = length(gamma);
    gamma_sol = gamma;
    srapp_prev = 0;
    ln_prev = 0;
    srapp_nxt = SRgapprox_active(p, G, H, gamma, gamma, sigma_sq, sigma_sqris, K, NR, BW);
    ln_nxt = srapp_nxt / func_gtot_active(p, gamma, H, sigma_sqris, K, Pca, mu);
    cvx_status = 'Failed';
    status = 'Solved';


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
    tol = 1e-8;
    
    while ln_nxt - ln_prev > tol || isempty(regexpi(cvx_status, 'Solved', 'once'))

        ng_iter = ng_iter + 1;
        srapp_prev = srapp_nxt;
        ln_prev = ln_nxt;
        gamma_prev = gamma_sol;

            
        cvx_begin
            cvx_quiet(true)
            cvx_precision high
            cvx_solver mosek
            variable gamma_opt(N) complex; % semidefinite % or hermitian semidefinite
            
            % GEE Dinkelbach objective: concave numerator - lambda * convex denominator.
            % Two points keep this a valid concave maximization:
            %   (1) the denominator uses func_gtot_active_cvx (quad_form), so CVX
            %       certifies the quadratic-in-gamma term as convex;
            %   (2) lambda is clamped to >= 0 with max(ln_prev,0): the rate
            %       surrogate (and hence the GEE estimate ln_prev) can be
            %       transiently negative for a poor intermediate gamma, and a
            %       negative coefficient would flip the convex denominator to
            %       concave, yielding an illegal "concave - concave" objective.
            lambda = opt_bool * max(ln_prev, 0);
            objective = SRgapprox_active_cvx(p, G, H, gamma_opt, gamma_prev, sigma_sq, sigma_sqris, K, NR, BW) - lambda*func_gtot_active_cvx(p, gamma_opt, H, sigma_sqris, K, Pca, mu);

            maximize objective;
            subject to
                if strcmpi(cons_mode, 'local')
                    abs(gamma_opt) <= 1;             % per-element |gamma_n| <= 1
                else
                    real(gamma_prev'*Rnorm*gamma_prev) + 2*real(gamma_prev'*Rnorm*(gamma_opt-gamma_prev)) >= trace(Rnorm);
                    quad_form(gamma_opt, Rnorm) <= PRmax;   % convex global budget (Eq. 7)
                end

                
        cvx_end

        % Robustness guard: if the solver did not return a usable (finite,
        % "Solved") point, keep the last good iterate and stop, rather than
        % propagating NaN/Inf into the next MMSE filter and CVX expression.
        if isempty(regexpi(cvx_status, 'Solved', 'once')) || any(~isfinite(gamma_opt))
            warning('gamma_cvxopt_active:solveFallback', 'CVX returned "%s" at iteration %d; keeping previous iterate.', cvx_status, ng_iter);
            gamma_sol = gamma_prev;
            cvx_status = 'Solved';   % let the while-loop terminate cleanly
            break;
        end

        gamma_sol = gamma_opt;
        srapp_nxt = SRgapprox_active(p, G, H, gamma_sol, gamma_prev, sigma_sq, sigma_sqris, K, NR, BW);
        ln_nxt = srapp_nxt / func_gtot_active(p, gamma_sol, H, sigma_sqris, K, Pca, mu);


        if opt_bool == 0
            ln_nxt =  srapp_nxt;    
        end

%         gamma_norm = norm(gamma_sol)^2;
%         fprintf('Step: %d, srapp_nxt: %f, ln_nxt: %f, gamma_norm: %f, Status of Solver: %s\n', ng_iter, srapp_nxt, ln_nxt, gamma_norm, cvx_status)
   end   
    
    srapp_sol = srapp_nxt;
    ln_sol = ln_nxt;

    if ln_nxt < ln_prev
        fprintf('True, ln_nxt is less than ln_prev \n')
        gamma_sol = gamma_prev;
        srapp_sol = srapp_prev;
        ln_sol = ln_prev;
    end
    
     if opt_bool == 0
        ln_sol = srapp_sol / func_gtot_active(p, gamma_sol, H, sigma_sqris, K, Pca, mu);
     end

    gamma_norm = norm(gamma_sol)^2;
    fprintf('1st approach -> Active Gamma Opt: Number of Steps: %d, SR_opt: %f, GEE_opt: %f, Gamma_norm: %f, Status of Solver: %s\n',ng_iter, srapp_sol, ln_sol, gamma_norm, cvx_status)
end
