function [p_opt, gamma_opt, ur_opt, sr_opt, gee_opt, npg_iter] = alt_opt1(p, pmax, G, H, gamma, sigma_sq, K, NR, PR, BW, mu, Pc, opt_bool, cons_mode)
% ALT_OPT1  ALGORITHM 3 (nearly-passive): alternating optimization of p, gamma and MMSE filters C.
%
%   function [p_opt, gamma_opt, ur_opt, sr_opt, gee_opt, npg_iter] = alt_opt1(p, pmax, G, H, gamma, sigma_sq, K, NR, PR, BW, mu, Pc, opt_bool, cons_mode)
%
% cons_mode (optional): 'global' (default) or 'local', forwarded to the RIS
%   subproblem gamma_cvxopt to select the reflection constraint (see Fig. 6).
%
% Paper reference: Algorithm 3, Problem (14)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    if nargin < 14 || isempty(cons_mode); cons_mode = 'global'; end

    mu0 = mu;
    Pc0 = Pc;
    
    if opt_bool == 0
        mu = 0;
        Pc = 1;
    end

    % Initialization
    gee_prev = 0;
    p_nxt = p; % initializing optimal power
    gamma_nxt = gamma; % initializing gamma_opt
    npg_iter = 0; % iteration

    % Linear MMSE update
    C_init = LMMSE_receiver_passive(p, G, H, gamma, sigma_sq, K, NR);

    % computing the next user / sum rate
    sinr_init = SINR_passive(p, C_init, G, H, gamma, sigma_sq, K);
    [ur_nxt, sr_nxt] = data_rate(sinr_init, K, BW);
    gee_nxt = sr_nxt / func_ptot(p, K, Pc, mu);
    
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

         % optimizing p
        [p_nxt, ~, ~, ~] = p_cvxopt(p_prev, pmax, G, H, gamma_nxt, sigma_sq, K, NR, BW, mu, Pc, opt_bool);
%         [p_opt, ~, ~] = optimization_p(p_opt, pmax, c_update, G, H, gamma_opt, sigma_sq, K, BW);
        
        
        % optimizing gamma
        [gamma_nxt, ~, ~, ~] = gamma_cvxopt(p_nxt, G, H, gamma_prev, sigma_sq, K, NR, PR, BW, mu, Pc, cons_mode);
        % [gamma_opt, ~, ~, ~] = optimization_gamma(p_opt, c_update, G, H, gamma_opt, sigma_sq, K, BW);

        % LMMSE Filter Optimization
        C_nxt = LMMSE_receiver_passive(p_nxt, G, H, gamma_nxt, sigma_sq, K, NR);
     
        % computing the actual sinr, user rate and sum rate using the
        % optimized values of (p, gamma) ~ (p_opt, gamma_opt)
        sinr_nxt = SINR_passive(p_nxt, C_nxt, G, H, gamma_nxt, sigma_sq, K);
        [ur_nxt, sr_nxt] = data_rate(sinr_nxt, K, BW);
        gee_nxt = sr_nxt / func_ptot(p_nxt, K, Pc, mu);   
        
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
        gee_opt = sr_opt / func_ptot(p_opt, K, Pc0, mu0);
    end

    fprintf('1st Approach -> Passive Alternating Opt: Number of Steps: %d, SR_opt: %f,  GEE_opt: %f\n', npg_iter, sr_opt, gee_opt)
end
