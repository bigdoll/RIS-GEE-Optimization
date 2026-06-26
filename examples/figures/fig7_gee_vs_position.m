% FIG7_GEE_VS_POSITION  Reproduce Fig. 7: achieved GEE (Algorithm 3) versus the
%                       RIS-to-BS distance d_{RIS-BS}, for two path-loss regimes.
%
% The end-to-end channel is the product of the RIS-BS and user-RIS links, so the
% optimal RIS placement depends on which link is more attenuated:
%
%   Case 1 (n_h = 4 > n_g = 2): the user-RIS link decays faster, so the RIS
%           should be placed closer to the users (larger d_{RIS-BS}).
%   Case 2 (n_h = 2 < n_g = 4): the RIS-BS link decays faster, so the RIS
%           should be placed closer to the BS (smaller d_{RIS-BS}).
%
% The per-link exponents are passed to generate_channels(..., nh, ng).
%
% Requirements: MATLAB + CVX + MOSEK.

clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(genpath(fullfile(thisDir, '..', '..', 'src')));

%% Parameters (Sec. V) --------------------------------------------------
K = 4; NR = 4; N = 100; BW0 = 20; mu = 1;
APPROACH = 1;                  % Algorithm 3
PR_active = 2;
Pmax = 10^(10/10);             % fixed transmit budget: 10 dBW
Pi   = (Pmax/K) * ones(K,1);
sigma_sq    = noise_power(-174, 10, BW0*1e6);
sigma_sqris = noise_power(-174, 10, BW0*1e6);
Pca = static_power_consumption(40, 20, 30, N);

r_cell = 100; hbs = 10; hris = 15; hmax_ue = 5;
posix_list = 10:15:85;         % RIS ground distance from the BS [m]
Ncarlo = 5;

% Two path-loss regimes: [nh (user-RIS), ng (RIS-BS)]
cases = struct('name', {'n_h=4, n_g=2 (RIS near users)', 'n_h=2, n_g=4 (RIS near BS)'}, ...
               'nh',   {4, 2}, ...
               'ng',   {2, 4});

%% Monte Carlo sweep ----------------------------------------------------
gee = zeros(numel(posix_list), numel(cases));

for ci = 1:numel(cases)
    for j = 1:numel(posix_list)
        ris_posix = posix_list(j);
        for mc = 1:Ncarlo
            [G, H] = generate_channels(r_cell, hbs, hris, hmax_ue, ris_posix, ...
                                       N, K, NR, cases(ci).nh, cases(ci).ng);
            phi = exp(1i*2*pi*rand(N,1));
            gee(j, ci) = gee(j, ci) + ...
                gee_active(APPROACH, Pi, Pmax, G, H, sigma_sq, sigma_sqris, ...
                           K, NR, PR_active, BW0, mu, Pca, phi);
        end
    end
end
gee = gee / Ncarlo;

%% Plot -----------------------------------------------------------------
figure;
plot(posix_list, gee(:,1), '-or', posix_list, gee(:,2), '-sm', 'LineWidth', 2);
grid on;
xlabel('RIS-to-BS distance  d_{RIS-BS} [m]');
ylabel('Global energy efficiency [Mbit/J]');
legend(cases(1).name, cases(2).name, 'Location', 'best');
title('Fig. 7 - GEE vs. RIS position (two path-loss regimes)');
