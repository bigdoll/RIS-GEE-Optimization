% RUN_ACTIVE_RIS  Reproduce the GEE / sum-rate vs. transmit-power curves for an
%                 ACTIVE RIS with a global reflection constraint.
%
% This script sweeps the maximum available transmit power and, for each value,
% evaluates five resource-allocation schemes (matching the paper's figures):
%
%   (1) Uniform power + random RIS allocation        (baseline, no optimization)
%   (2) Approach 1, sum-rate optimal   (Algorithm 3, opt_bool = 0)
%   (3) Approach 1, GEE optimal        (Algorithm 3, opt_bool = 1)
%   (4) Approach 2, sum-rate optimal   (Algorithm 6, opt_bool = 0)
%   (5) Approach 2, GEE optimal        (Algorithm 6, opt_bool = 1)
%
% Approach 1 = alternating optimization of (p, gamma, C)          [Sec. IV-A]
% Approach 2 = MMSE filters embedded, optimization of (p, X), SDR [Sec. IV-B]
%
% A "warm-start" initialization is used: the optimum found for one power level
% initializes the next level, which speeds up convergence (paper, Sec. V).
%
% Requirements: MATLAB + CVX (http://cvxr.com/cvx) with the MOSEK solver.
%
% Paper: R. K. Fotock, A. Zappone, M. Di Renzo, "Energy Efficiency Optimization
% in RIS-Aided Wireless Networks: Active Versus Nearly-Passive RIS With Global
% Reflection Constraints," IEEE Trans. Commun., 2024.

clear; clc; close all;

% Put every function of the repository on the MATLAB path.
thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(thisDir, '..', 'src')));

% Reset any CVX state left over from a previously interrupted run, so the
% first cvx_begin does not warn "a non-empty cvx problem already exists".
if exist('cvx_clear', 'file'); cvx_clear; end

%% ------------------------------------------------------------------ %%
%%  1. SYSTEM PARAMETERS  (values from Sec. V of the paper)            %%
%% ------------------------------------------------------------------ %%
K   = 4;        % number of single-antenna users
NR  = 4;        % number of BS receive antennas
N   = 100;      % number of RIS reflecting elements
BW  = 20e6;     % communication bandwidth [Hz]
BW0 = BW / 1e6; % bandwidth in MHz (rates are reported in Mbit/s)
mu  = 1;        % inverse efficiency of the users' transmit amplifiers (mu >= 1)

% RIS regime:  PR = 1 -> nearly-passive,  PR > 1 -> active amplification budget.
PR = 2;

% Transmit-power sweep  Pmax = P_step grid from p_min to p_max [dBW]
p_min = -40; p_max = 40; p_step = 2;
[P_Watts, Ni] = available_power(p_min, p_max, p_step);

% Static power consumption of an ACTIVE RIS [dBm]  (Sec. V)
Po_dBm     = 40;   % all nodes except the RIS                (P0)
Pcn_dBm    = 20;   % per active RIS element                  (Pac,n)
Po_ris_dBm = 30;   % other static sources at the active RIS  (P0,RIS^(a))
Pca = static_power_consumption(Po_dBm, Pcn_dBm, Po_ris_dBm, N);

% Receiver and RIS noise   (PSD -174 dBm/Hz, noise figure 10 dB)
noise_psd = -174; noise_figure = 10;
sigma_sq    = noise_power(noise_psd, noise_figure, BW);   % receiver noise variance
sigma_sqris = noise_power(noise_psd, noise_figure, BW);   % RIS amplification noise

% Network geometry  (Sec. V)
r_cell    = 100;   % users scattered within 100 m of the RIS
hbs       = 10;    % BS height  [m]
hris      = 15;    % RIS height [m]
hmax_ue   = 5;     % maximum user height [m]
ris_posix = 50;    % BS-to-RIS ground distance [m]

% Monte Carlo channel realizations
Ncarlo = 10;

%% ------------------------------------------------------------------ %%
%%  2. CHANNEL GENERATION                                              %%
%% ------------------------------------------------------------------ %%
% G: NR x N  RIS -> BS channels ;  H: N x K  users -> RIS channels.
G = zeros(NR, N, Ncarlo);
H = zeros(N,  K, Ncarlo);
for mc = 1:Ncarlo
    [G(:,:,mc), H(:,:,mc)] = generate_channels(r_cell, hbs, hris, hmax_ue, ris_posix, N, K, NR);
end
% NOTE: this example uses the first realization (mc = 1). Wrap the loop below
% in an outer "for mc = 1:Ncarlo" and average to reproduce the paper's curves.
G = G(:,:,1);  H = H(:,:,1);

%% ------------------------------------------------------------------ %%
%%  3. PRE-ALLOCATION                                                  %%
%% ------------------------------------------------------------------ %%
SR1      = zeros(Ni,1);  GEE1      = zeros(Ni,1);   % uniform baseline
SR1_ropt = zeros(Ni,1);  GEE1_ropt = zeros(Ni,1);   % approach 1, sum-rate opt
SR1_eopt = zeros(Ni,1);  GEE1_eopt = zeros(Ni,1);   % approach 1, GEE opt
SR2_ropt = zeros(Ni,1);  GEE2_ropt = zeros(Ni,1);   % approach 2, sum-rate opt
SR2_eopt = zeros(Ni,1);  GEE2_eopt = zeros(Ni,1);   % approach 2, GEE opt
gamma_ropt = zeros(N,Ni); gamma_eopt = zeros(N,Ni);
X_ropt = zeros(N,N,Ni);   X_eopt = zeros(N,N,Ni);

opt_SR  = 0;   % opt_bool = 0 -> maximize sum-rate
opt_GEE = 1;   % opt_bool = 1 -> maximize GEE

phi = exp(1i*2*pi*rand(N,1));   % random RIS phases (shared feasible start)

%% ------------------------------------------------------------------ %%
%%  4. POWER SWEEP                                                     %%
%% ------------------------------------------------------------------ %%
for i = 1:Ni

    Pi_max = P_Watts(i);
    Pi     = (Pi_max / K) * ones(K,1);   % uniform power allocation
    fprintf('\n=== Pmax = %g W (point %d/%d) ===\n', Pi_max, i, Ni);

    % Feasible RIS start: scale the random phases to satisfy the global
    % reflection budget  tr(R) <= tr(R*gamma*gamma') <= PRmax + tr(R).
    R     = func_R(Pi, H, sigma_sqris, K, N);
    Rnorm = R / abs(max(R, [], 'all'));
    PRmax = sqrt((1/4) * trace(Rnorm) * N * PR);
    rho   = sqrt(PRmax / abs(phi' * Rnorm * phi));
    gamma = rho .* phi;
    X     = gamma * gamma';

    % Warm start: reuse the previous power level's optimum.
    if i == 1
        P1_rin = Pi; P1_ein = Pi; P2_rin = Pi; P2_ein = Pi;
        gamma_rin = gamma; gamma_ein = gamma; X_rin = X; X_ein = X;
    end

    % (1) Uniform baseline
    SR1(i)  = SRp_active(Pi, G, H, gamma, sigma_sq, sigma_sqris, K, NR, BW0);
    GEE1(i) = SR1(i) / func_ptot_active(Pi, gamma, H, sigma_sqris, K, Pca, mu);

    % (2)-(3) Approach 1 - Algorithm 3
    [P1_rin, gamma_ropt(:,i), ~, SR1_ropt(i), GEE1_ropt(i)] = ...
        altopt1_active(P1_rin, Pi_max, G, H, gamma_rin, sigma_sq, sigma_sqris, K, NR, PR, BW0, mu, Pca, opt_SR);
    gamma_rin = gamma_ropt(:,i);

    [P1_ein, gamma_eopt(:,i), ~, SR1_eopt(i), GEE1_eopt(i)] = ...
        altopt1_active(P1_ein, Pi_max, G, H, gamma_ein, sigma_sq, sigma_sqris, K, NR, PR, BW0, mu, Pca, opt_GEE);
    gamma_ein = gamma_eopt(:,i);

    % (4)-(5) Approach 2 - Algorithm 6
    [P2_rin, X_ropt(:,:,i), SR2_ropt(i), GEE2_ropt(i)] = ...
        altopt2_active(P2_rin, Pi_max, G, H, X_rin, sigma_sq, sigma_sqris, K, NR, PR, BW0, mu, Pca, opt_SR);
    X_rin = X_ropt(:,:,i);

    [P2_ein, X_eopt(:,:,i), SR2_eopt(i), GEE2_eopt(i)] = ...
        altopt2_active(P2_ein, Pi_max, G, H, X_ein, sigma_sq, sigma_sqris, K, NR, PR, BW0, mu, Pca, opt_GEE);
    X_ein = X_eopt(:,:,i);
end

%% ------------------------------------------------------------------ %%
%%  5. PLOTS                                                           %%
%% ------------------------------------------------------------------ %%
P_dBW = p_min:p_step:p_max;

figure;
plot(P_dBW, SR1,'-ob', P_dBW, SR1_ropt,'-dk', P_dBW, SR2_ropt,'-og', ...
     P_dBW, SR1_eopt,'-ro', P_dBW, SR2_eopt,'-mo', 'LineWidth', 2);
grid on;
legend('Uniform (no opt.)', 'Approach 1 - SR opt', 'Approach 2 - SR opt', ...
       'Approach 1 - GEE opt', 'Approach 2 - GEE opt', 'Location', 'best');
xlabel('Maximum available power P_{max} [dBW]');
ylabel('Sum rate [Mbit/s]');
title('Active RIS: sum rate vs. transmit power');

figure;
plot(P_dBW, GEE1,'-ob', P_dBW, GEE1_ropt,'-dk', P_dBW, GEE2_ropt,'-og', ...
     P_dBW, GEE1_eopt,'-ro', P_dBW, GEE2_eopt,'-mo', 'LineWidth', 2);
grid on;
legend('Uniform (no opt.)', 'Approach 1 - SR opt', 'Approach 2 - SR opt', ...
       'Approach 1 - GEE opt', 'Approach 2 - GEE opt', 'Location', 'best');
xlabel('Maximum available power P_{max} [dBW]');
ylabel('Global energy efficiency [Mbit/J]');
title('Active RIS: GEE vs. transmit power');
