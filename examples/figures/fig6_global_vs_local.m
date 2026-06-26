% FIG6_GLOBAL_VS_LOCAL  Reproduce Fig. 6: achieved GEE versus transmit power,
%                       comparing the GLOBAL reflection constraint with the
%                       classical LOCAL (per-element |gamma_n| <= 1) constraint.
%
% Both curves use Algorithm 3 (alternating optimization). The only difference
% is the RIS feasibility set, selected through the cons_mode argument that the
% drivers/solvers now accept ('global' vs 'local'). The paper shows that the
% global constraint always performs at least as well as the local one.
%
% Requirements: MATLAB + CVX + MOSEK.

clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(genpath(fullfile(thisDir, '..', '..', 'src')));

%% Parameters (Sec. V) --------------------------------------------------
K = 4; NR = 4; N = 100; BW0 = 20; mu = 1;
APPROACH = 1;                  % global-vs-local variant is implemented for Alg. 3
PR_active = 2; PR_passive = 1;
sigma_sq    = noise_power(-174, 10, BW0*1e6);
sigma_sqris = noise_power(-174, 10, BW0*1e6);
Pca = static_power_consumption(40, 20, 30, N);
Pc  = static_power_consumption(40,  0, 20, N);

p_dBW_list = -10:10:40;        % coarse transmit-power sweep [dBW]
Ncarlo = 5;
geo = {100, 10, 15, 5, 50};

%% Monte Carlo sweep ----------------------------------------------------
nP = numel(p_dBW_list);
[a_glob, a_loc, p_glob, p_loc] = deal(zeros(nP,1));

for mc = 1:Ncarlo
    [G, H] = generate_channels(geo{:}, N, K, NR);
    phi = exp(1i*2*pi*rand(N,1));
    for j = 1:nP
        Pmax = 10^(p_dBW_list(j)/10);
        Pi   = (Pmax/K) * ones(K,1);
        a_glob(j) = a_glob(j) + gee_active(APPROACH, Pi, Pmax, G, H, sigma_sq, sigma_sqris, K, NR, PR_active, BW0, mu, Pca, phi, 1, 'global');
        a_loc(j)  = a_loc(j)  + gee_active(APPROACH, Pi, Pmax, G, H, sigma_sq, sigma_sqris, K, NR, PR_active, BW0, mu, Pca, phi, 1, 'local');
        p_glob(j) = p_glob(j) + gee_passive(APPROACH, Pi, Pmax, G, H, sigma_sq, K, NR, PR_passive, BW0, mu, Pc, phi, 1, 'global');
        p_loc(j)  = p_loc(j)  + gee_passive(APPROACH, Pi, Pmax, G, H, sigma_sq, K, NR, PR_passive, BW0, mu, Pc, phi, 1, 'local');
    end
end
a_glob = a_glob/Ncarlo; a_loc = a_loc/Ncarlo;
p_glob = p_glob/Ncarlo; p_loc = p_loc/Ncarlo;

%% Plot -----------------------------------------------------------------
figure;
plot(p_dBW_list, a_glob, '-or', p_dBW_list, a_loc, '--xr', ...
     p_dBW_list, p_glob, '-sb', p_dBW_list, p_loc, '--db', 'LineWidth', 2);
grid on;
xlabel('Maximum available power P_{max} [dBW]');
ylabel('Global energy efficiency [Mbit/J]');
legend('Active - global', 'Active - local', ...
       'Passive - global', 'Passive - local', 'Location', 'best');
title('Fig. 6 - global vs. local reflection constraint');
