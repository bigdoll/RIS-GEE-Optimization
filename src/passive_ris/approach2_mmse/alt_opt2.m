function [p_sol, X_sol, sr_sol, gee_sol, Npg_iter] = alt_opt2(p, pmax, G, H, X, sigma_sq, K, NR, PR, BW, mu, Pc, opt_bool)
% ALT_OPT2  ALGORITHM 6 (nearly-passive): alternating optimization of X and p (MMSE embedded).
%
%   function [p_sol, X_sol, sr_sol, gee_sol, Npg_iter] = alt_opt2(p, pmax, G, H, X, sigma_sq, K, NR, PR, BW, mu, Pc, opt_bool)
%
% Paper reference: Algorithm 6, Problem (47)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    mu0 = mu;
    Pc0 = Pc;
    
    if opt_bool == 0
        mu = 0;
        Pc = 1;
    end

    % Initialization
    sr_prev = 0; % initializing previous sum rate
    gee_prev = 0; % initializing previous GEE
    p_opt = p; % initializing optimal power
    X_opt = X; % initializing X_opt
    Npg_iter = 0; % iteration

    % computing the next sum rate and GEE
    sr_nxt = func_SRmmse(p, G, H, X, sigma_sq, K, NR, BW);
    gee_nxt = sr_nxt/func_ptot(p,K,Pc,mu);
    
    tol = 1e-6;% tolerance
    
    while gee_nxt - gee_prev > tol
        
        Npg_iter = Npg_iter + 1;
        
        sr_prev = sr_nxt;

        gee_prev = gee_nxt;

        p_prev = p_opt;

        X_prev = X_opt;
         
        % optimizing p
        [p_opt, ~, ~, ~] = power2_cvxopt(p_opt, pmax, G, H, X_opt, sigma_sq, K, NR, BW, mu, Pc, opt_bool);
       
        % optimizing X
        [X_opt, ~, ~, ~] = cvxopt_X(p_opt, G, H, X_opt, sigma_sq, K, NR, PR, BW, mu, Pc);
         
        % computing next SR and GEE
        sr_nxt = func_SRmmse(p_opt, G, H, X_opt, sigma_sq, K, NR, BW);
        gee_nxt = sr_nxt / func_ptot(p_opt,K,Pc,mu);
        
    end
    
    p_sol = p_opt;
    X_sol = X_opt;
    sr_sol = sr_nxt;
    gee_sol = gee_nxt;

    if gee_nxt < gee_prev
        p_sol = p_prev;
        X_sol = X_prev;
        sr_sol = sr_prev;
        gee_sol = gee_prev;
    end
    
    if opt_bool == 0
        gee_sol = sr_sol / func_ptot(p_sol, K, Pc0, mu0);
    end
 
    fprintf('2nd Approach -> Passive Alternating Opt: Number of Steps: %d, SR_opt: %f, GEE_opt: %f\n', Npg_iter, sr_sol, gee_sol);
end
