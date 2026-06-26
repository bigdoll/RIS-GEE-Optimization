function [p_sol, X_sol, sr_sol, gee_sol, Npg_iter] = altopt2_active(p, pmax, G, H, X, sigma_sq, sigma_sqris, K, NR, PR, BW, mu, Pca, opt_bool)
% ALTOPT2_ACTIVE  ALGORITHM 6 (active): alternating optimization of X = gamma*gamma^H and p (MMSE embedded).
%
%   function [p_sol, X_sol, sr_sol, gee_sol, Npg_iter] = altopt2_active(p, pmax, G, H, X, sigma_sq, sigma_sqris, K, NR, PR, BW, mu, Pca, opt_bool)
%
% Paper reference: Algorithm 6, Problem (47)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    % Initialization
    sr_prev = 0; % initializing previous sum rate
    gee_prev = 0; % initializing previous GEE
    p_opt = p; % initializing optimal power
    X_opt = X; % initializing X_opt
    Npg_iter = 0; % iteration

    % computing the next sum rate and GEE
    sr_nxt = SRX_active2(p, G, H, X, sigma_sq, sigma_sqris, K, NR, BW);
    gee_nxt = sr_nxt/func_Xtot_active2(p, X, H, sigma_sqris, K, Pca, mu);
 
     if opt_bool == 0
       gee_nxt =  sr_nxt;    
     end
    
    tol = 1e-6;% tolerance
    
    while gee_nxt - gee_prev > tol
        
        Npg_iter = Npg_iter + 1;
        
        sr_prev = sr_nxt;

        gee_prev = gee_nxt;

        p_prev = p_opt;

        X_prev = X_opt;

        % optimizing X
        [X_opt, ~, ~, ~] = X_cvxopt_active2(p_opt, G, H, X_opt, sigma_sq, sigma_sqris, K, NR, PR, BW, mu, Pca, opt_bool); 

         % optimizing p
        [p_opt, ~, ~, ~] = p_cvxopt_active2(p_opt, pmax, G, H, X_opt, sigma_sq, sigma_sqris, K, NR, BW, mu, Pca, opt_bool);
               
        % computing next SR and GEE
        sr_nxt = SRX_active2(p_opt, G, H, X_opt, sigma_sq, sigma_sqris, K, NR, BW);
        gee_nxt = sr_nxt / func_Xtot_active2(p_opt, X_opt, H, sigma_sqris, K, Pca, mu);
        
         if opt_bool == 0
            gee_nxt =  sr_nxt;    
         end
        
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
        gee_sol = sr_sol / func_Xtot_active2(p_sol, X_sol, H, sigma_sqris, K, Pca, mu);
    end
    
    fprintf('2nd Approach -> Alternating Opt: Number of Steps: %d, SR_opt: %f, GEE_opt: %f\n', Npg_iter, sr_sol, gee_sol);
end
