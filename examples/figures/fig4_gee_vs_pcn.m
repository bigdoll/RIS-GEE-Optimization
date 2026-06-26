% FIG4_GEE_VS_PCN  Reproduce Fig. 4: achieved GEE versus the per-element static
%                  power of an ACTIVE RIS, Pc,n^(a), compared with a
%                  nearly-passive RIS (whose GEE does not depend on Pc,n^(a)).
%
% This sweep reveals the operating regime in which an active RIS is more
% energy-efficient than a nearly-passive one: as the active hardware power
% Pc,n^(a) grows, the active-RIS GEE falls below the passive curve.
%
% Requirements: MATLAB + CVX + MOSEK.   Runtime grows with Ncarlo and the sweep.

clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);                                   % gee_active / gee_passive
addpath(genpath(fullfile(thisDir, '..', '..', 'src')));

%% Parameters (Sec. V) --------------------------------------------------
K = 4; NR = 4; N = 100; BW0 = 20; mu = 1;
APPROACH = 1;                  % 1 = Algorithm 3, 2 = Algorithm 6
PR_active = 2; PR_passive = 1;
Pmax  = 10^(10/10);            % fixed transmit budget: 10 dBW
Pi    = (Pmax/K) * ones(K,1);
sigma_sq    = noise_power(-174, 10, BW0*1e6);
sigma_sqris = noise_power(-174, 10, BW0*1e6);

pcn_dBm_list = 0:5:40;         % active per-element static power sweep [dBm]
Ncarlo = 5;                    % increase for smoother curves

% Geometry
geo = {100, 10, 15, 5, 50};    % r_cell, hbs, hris, hmax_ue, ris_posix

%% Monte Carlo sweep ----------------------------------------------------
gee_active_v  = zeros(numel(pcn_dBm_list), 1);
gee_passive_v = zeros(numel(pcn_dBm_list), 1);

for mc = 1:Ncarlo
    [G, H] = generate_channels(geo{:}, N, K, NR);
    phi = exp(1i*2*pi*rand(N,1));

    % Nearly-passive RIS: static power independent of Pc,n^(a).
    Pc = static_power_consumption(40, 0, 20, N);
    gee_passive_v = gee_passive_v + ...
        gee_passive(APPROACH, Pi, Pmax, G, H, sigma_sq, K, NR, PR_passive, BW0, mu, Pc, phi);

    % Active RIS: recompute static power for each Pc,n^(a).
    for j = 1:numel(pcn_dBm_list)
        Pca = static_power_consumption(40, pcn_dBm_list(j), 30, N);
        gee_active_v(j) = gee_active_v(j) + ...
            gee_active(APPROACH, Pi, Pmax, G, H, sigma_sq, sigma_sqris, K, NR, PR_active, BW0, mu, Pca, phi);
    end
end
gee_active_v  = gee_active_v  / Ncarlo;
gee_passive_v = (gee_passive_v / Ncarlo) * ones(size(pcn_dBm_list))';  % flat line

%% Plot -----------------------------------------------------------------
figure;
plot(pcn_dBm_list, gee_active_v, '-or', pcn_dBm_list, gee_passive_v, '-sb', 'LineWidth', 2);
grid on;
xlabel('Active RIS per-element static power  P_{c,n}^{(a)} [dBm]');
ylabel('Global energy efficiency [Mbit/J]');
legend('Active RIS', 'Nearly-passive RIS', 'Location', 'best');
title('Fig. 4 - GEE vs. active RIS static power');
