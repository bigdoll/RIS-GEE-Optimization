% FIG5_GEE_VS_N  Reproduce Fig. 5: achieved GEE versus the number of RIS
%                reflecting elements N, for an active and a nearly-passive RIS.
%
% N is taken equal for both RIS types; the active per-element static power is
% fixed to Pc,n^(a) = 20 dBm (paper caption). Channels and static powers are
% regenerated for every value of N.
%
% Requirements: MATLAB + CVX + MOSEK.

clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(genpath(fullfile(thisDir, '..', '..', 'src')));

%% Parameters (Sec. V) --------------------------------------------------
K = 4; NR = 4; BW0 = 20; mu = 1;
APPROACH = 1;                  % 1 = Algorithm 3, 2 = Algorithm 6
PR_active = 2; PR_passive = 1;
Pmax = 10^(10/10);             % fixed transmit budget: 10 dBW
Pi   = (Pmax/K) * ones(K,1);
sigma_sq    = noise_power(-174, 10, BW0*1e6);
sigma_sqris = noise_power(-174, 10, BW0*1e6);

N_list = [50 100 150 200];     % RIS size sweep
Ncarlo = 5;
geo = {100, 10, 15, 5, 50};

%% Monte Carlo sweep ----------------------------------------------------
gee_active_v  = zeros(numel(N_list), 1);
gee_passive_v = zeros(numel(N_list), 1);

for j = 1:numel(N_list)
    N = N_list(j);
    Pca = static_power_consumption(40, 20, 30, N);   % active,  Pc,n^(a) = 20 dBm
    Pc  = static_power_consumption(40,  0, 20, N);   % nearly-passive

    for mc = 1:Ncarlo
        [G, H] = generate_channels(geo{:}, N, K, NR);
        phi = exp(1i*2*pi*rand(N,1));
        gee_active_v(j)  = gee_active_v(j)  + ...
            gee_active(APPROACH, Pi, Pmax, G, H, sigma_sq, sigma_sqris, K, NR, PR_active, BW0, mu, Pca, phi);
        gee_passive_v(j) = gee_passive_v(j) + ...
            gee_passive(APPROACH, Pi, Pmax, G, H, sigma_sq, K, NR, PR_passive, BW0, mu, Pc, phi);
    end
end
gee_active_v  = gee_active_v  / Ncarlo;
gee_passive_v = gee_passive_v / Ncarlo;

%% Plot -----------------------------------------------------------------
figure;
plot(N_list, gee_active_v, '-or', N_list, gee_passive_v, '-sb', 'LineWidth', 2);
grid on;
xlabel('Number of RIS reflecting elements  N');
ylabel('Global energy efficiency [Mbit/J]');
legend('Active RIS', 'Nearly-passive RIS', 'Location', 'best');
title('Fig. 5 - GEE vs. number of RIS elements');
