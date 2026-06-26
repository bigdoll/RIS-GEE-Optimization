% TABLES_CONVERGENCE  Reproduce Tables I and II: average convergence time and
%                     average number of iterations for Algorithms 3 and 6, with
%                     an active and a nearly-passive RIS.
%
% Results are averaged over Ncarlo channel realizations and a set of transmit
% powers, mirroring the paper's methodology (convergence tolerance ~ 1e-3..1e-6,
% set inside the drivers). Wall-clock times depend on your machine and on the
% CVX/MOSEK installation, so absolute values will differ from the paper's
% workstation; the relative trends (Alg. 6 slower than Alg. 3; active slower
% than passive) should hold.
%
% Requirements: MATLAB + CVX + MOSEK.

clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(genpath(fullfile(thisDir, '..', '..', 'src')));

%% Parameters (Sec. V) --------------------------------------------------
K = 4; NR = 4; N = 100; BW0 = 20; mu = 1;
PR_active = 2; PR_passive = 1;
sigma_sq    = noise_power(-174, 10, BW0*1e6);
sigma_sqris = noise_power(-174, 10, BW0*1e6);
Pca = static_power_consumption(40, 20, 30, N);
Pc  = static_power_consumption(40,  0, 20, N);
geo = {100, 10, 15, 5, 50};

p_dBW_list = [0 20];           % transmit powers to average over [dBW]
Ncarlo = 5;

%% Accumulators  [Alg3_passive Alg6_passive Alg3_active Alg6_active] -----
t_sum    = zeros(1,4);
iter_sum = zeros(1,4);
nRuns    = 0;

for mc = 1:Ncarlo
    [G, H] = generate_channels(geo{:}, N, K, NR);
    phi = exp(1i*2*pi*rand(N,1));
    for pp = p_dBW_list
        Pmax = 10^(pp/10);
        Pi   = (Pmax/K) * ones(K,1);

        [~,~,it,t] = gee_passive(1, Pi, Pmax, G, H, sigma_sq, K, NR, PR_passive, BW0, mu, Pc, phi);
        t_sum(1)=t_sum(1)+t; iter_sum(1)=iter_sum(1)+it;
        [~,~,it,t] = gee_passive(2, Pi, Pmax, G, H, sigma_sq, K, NR, PR_passive, BW0, mu, Pc, phi);
        t_sum(2)=t_sum(2)+t; iter_sum(2)=iter_sum(2)+it;
        [~,~,it,t] = gee_active(1, Pi, Pmax, G, H, sigma_sq, sigma_sqris, K, NR, PR_active, BW0, mu, Pca, phi);
        t_sum(3)=t_sum(3)+t; iter_sum(3)=iter_sum(3)+it;
        [~,~,it,t] = gee_active(2, Pi, Pmax, G, H, sigma_sq, sigma_sqris, K, NR, PR_active, BW0, mu, Pca, phi);
        t_sum(4)=t_sum(4)+t; iter_sum(4)=iter_sum(4)+it;

        nRuns = nRuns + 1;
    end
end
t_avg    = t_sum    / nRuns;
iter_avg = iter_sum / nRuns;

%% Print ----------------------------------------------------------------
fprintf('\n==================  TABLE I  -  average convergence time [s]  ==================\n');
fprintf('%-18s %12s %12s\n', '', 'Algorithm 3', 'Algorithm 6');
fprintf('%-18s %12.3f %12.3f\n', 'Nearly-passive', t_avg(1), t_avg(2));
fprintf('%-18s %12.3f %12.3f\n', 'Active',         t_avg(3), t_avg(4));

fprintf('\n==================  TABLE II -  average number of iterations  ==================\n');
fprintf('%-18s %12s %12s\n', '', 'Algorithm 3', 'Algorithm 6');
fprintf('%-18s %12.2f %12.2f\n', 'Nearly-passive', iter_avg(1), iter_avg(2));
fprintf('%-18s %12.2f %12.2f\n', 'Active',         iter_avg(3), iter_avg(4));
fprintf('\n(averaged over %d runs: %d realizations x %d power levels)\n', nRuns, Ncarlo, numel(p_dBW_list));
