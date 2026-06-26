function [gee, sr, iters, t] = gee_active(approach, Pi, Pi_max, G, H, sigma_sq, sigma_sqris, K, NR, PR, BW0, mu, Pca, phi, opt_bool, cons_mode)
% GEE_ACTIVE  Evaluate the GEE achieved by one active-RIS algorithm on one
%             channel realization. Helper shared by the figure scripts.
%
%   [gee, sr, iters, t] = gee_active(approach, Pi, Pi_max, G, H, sigma_sq, ...
%                            sigma_sqris, K, NR, PR, BW0, mu, Pca, phi, ...
%                            opt_bool, cons_mode)
%
%   approach  : 1 -> Algorithm 3 (alternating p, gamma, C)
%               2 -> Algorithm 6 (MMSE embedded, optimize p and X)
%   opt_bool  : 1 -> maximize GEE (default), 0 -> maximize sum-rate
%   cons_mode : 'global' (default) or 'local'  (approach 1 only; Fig. 6)
%   phi       : N x 1 unit-modulus phases used to build a feasible RIS start.
%
% Returns the achieved GEE, sum-rate, iteration count and wall-clock time.

    if nargin < 15 || isempty(opt_bool);  opt_bool  = 1;        end
    if nargin < 16 || isempty(cons_mode); cons_mode = 'global'; end

    N = numel(phi);

    % Feasible RIS start satisfying the global reflection budget.
    R     = func_R(Pi, H, sigma_sqris, K, N);
    Rnorm = R / abs(max(R, [], 'all'));
    PRmax = sqrt((1/4) * trace(Rnorm) * N * PR);
    rho   = sqrt(PRmax / abs(phi' * Rnorm * phi));
    gamma = rho .* phi;

    if approach == 1
        tic;
        [~, ~, ~, sr, gee, iters] = altopt1_active(Pi, Pi_max, G, H, gamma, ...
            sigma_sq, sigma_sqris, K, NR, PR, BW0, mu, Pca, opt_bool, cons_mode);
        t = toc;
    else
        X = gamma * gamma';
        tic;
        [~, ~, sr, gee, iters] = altopt2_active(Pi, Pi_max, G, H, X, ...
            sigma_sq, sigma_sqris, K, NR, PR, BW0, mu, Pca, opt_bool);
        t = toc;
    end
end
