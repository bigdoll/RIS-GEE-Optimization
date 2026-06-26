function [p_opt, gamma_opt, ur_opt, sr_opt, gee_opt, npg_iter] = altopt1_active(p, pmax, G, H, gamma, sigma_sq, sigma_sqris, K, NR, PR, BW, mu, Pca, opt_bool, cons_mode)
% ALTOPT1_ACTIVE  ALGORITHM 3 (active): alternating optimization of power p, RIS vector gamma and MMSE filters C.
%
%   function [p_opt, gamma_opt, ur_opt, sr_opt, gee_opt, npg_iter] = altopt1_active(p, pmax, G, H, gamma, sigma_sq, sigma_sqris, K, NR, PR, BW, mu, Pca, opt_bool, cons_mode)
%
% cons_mode (optional): 'global' (default) or 'local', forwarded to the RIS
%   subproblem gamma_cvxopt_active to select the reflection constraint (Fig. 6).
%
% Paper reference: Algorithm 3, Problem (12)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    if nargin < 15 || isempty(cons_mode); cons_mode = 'global'; end

    % Initialization
    gee_prev = 0;
    p_nxt = p; % initializing optimal power
    gamma_nxt = gamma; % initializing gamma_opt
    npg_iter = 0; % iteration

    % Linear MMSE update
    C_init = LMMSE_receiver_active(p, G, H, gamma, sigma_sq, sigma_sqris, K, NR);

    % computing the next user / sum rate
    sinr_init = SINR_active(p, C_init, G, H, gamma, sigma_sq, sigma_sqris, K);
    [ur_nxt, sr_nxt] = data_rate_active(sinr_init, K, BW);
    gee_nxt = sr_nxt / func_gtot_active(p, gamma, H, sigma_sqris, K, Pca, mu);
   
    if opt_bool == 0
       gee_nxt =  sr_nxt;    
    end

    tol = 1e-6;% tolerance
    
    while gee_nxt - gee_prev > tol 
        
        %count
        npg_iter = npg_iter + 1;
        
        % update previous values
        p_prev = p_nxt;
        gamma_prev = gamma_nxt;
        ur_prev = ur_nxt;
        sr_prev = sr_nxt; 
        gee_prev = gee_nxt;
        
        % optimizing gamma
        [gamma_nxt, ~, ~, ~] = gamma_cvxopt_active(p_nxt, G, H, gamma_prev, sigma_sq, sigma_sqris, K, NR, PR, BW, mu, Pca, opt_bool, cons_mode);
        % [gamma_opt, ~, ~, ~] = optimization_gamma(p_opt, c_update, G, H, gamma_opt, sigma_sq, K, BW);
        
         % optimizing p
        [p_nxt, ~, ~, ~] = p_cvxopt_active(p_prev, pmax, G, H, gamma_nxt, sigma_sq, sigma_sqris, K, NR, BW, mu, Pca, opt_bool);
%       p_opt, ~, ~] = optimization_p(p_opt, pmax, c_update, G, H, gamma_opt, sigma_sq, K, BW);
        
       
        % LMMSE Filter Optimization
        C_nxt = LMMSE_receiver_active(p_nxt, G, H, gamma_nxt, sigma_sq, sigma_sqris, K, NR);
     
        % computing the actual sinr, user rate and sum rate using the
        % optimized values of (p, gamma) ~ (p_opt, gamma_opt)
        sinr_nxt = SINR_active(p_nxt, C_nxt, G, H, gamma_nxt, sigma_sq, sigma_sqris, K);
        [ur_nxt, sr_nxt] = data_rate_active(sinr_nxt, K, BW);
        gee_nxt = sr_nxt / func_gtot_active(p_nxt, gamma_nxt, H, sigma_sqris, K, Pca, mu);   
        
         if opt_bool == 0
            gee_nxt =  sr_nxt;    
         end

    end
     
     p_opt = p_nxt;
     gamma_opt = gamma_nxt;
     ur_opt = ur_nxt;
     sr_opt = sr_nxt;
     gee_opt = gee_nxt;
    
    if gee_nxt < gee_prev
        p_opt = p_prev;
        gamma_opt = gamma_prev;
        ur_opt = ur_prev;
        sr_opt = sr_prev;
        gee_opt = gee_prev;
    end
    
    if opt_bool == 0
        gee_opt = sr_opt / func_gtot_active(p_opt, gamma_opt, H, sigma_sqris, K, Pca, mu);
    end

    fprintf('1st approach -> Active - Alternating Opt: Number of Steps: %d, SR_opt Value: %f, GEE_opt Value : %f\n', npg_iter, sr_opt, gee_opt)
end
