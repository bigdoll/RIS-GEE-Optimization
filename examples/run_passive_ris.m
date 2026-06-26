% RUN_PASSIVE_RIS  Reproduce the GEE / sum-rate vs. transmit-power curves for a
%                  NEARLY-PASSIVE RIS with a global reflection constraint.
%
% The nearly-passive RIS is the special case of the active RIS obtained by
% setting the amplification budget to zero and the RIS noise to zero
% (PR = 1, sigma_RIS = 0). See Problem (14) in the paper. The script evaluates
% the same five schemes as run_active_ris.m:
%
%   (1) Uniform power + random RIS allocation        (baseline, no optimization)
%   (2) Approach 1, sum-rate optimal   (Algorithm 3, opt_bool = 0)
%   (3) Approach 1, GEE optimal        (Algorithm 3, opt_bool = 1)
%   (4) Approach 2, sum-rate optimal   (Algorithm 6, opt_bool = 0)
%   (5) Approach 2, GEE optimal        (Algorithm 6, opt_bool = 1)
%
% Requirements: MATLAB + CVX (http://cvxr.com/cvx) with the MOSEK solver.
%
% Paper: R. K. Fotock, A. Zappone, M. Di Renzo, "Energy Efficiency Optimization
% in RIS-Aided Wireless Networks: Active Versus Nearly-Passive RIS With Global
% Reflection Constraints," IEEE Trans. Commun., 2024.

clear; clc; close all;

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
BW0 = BW / 1e6; % bandwidth in MHz (rates reported in Mbit/s)
mu  = 1;        % inverse efficiency of the users' transmit amplifiers (mu >= 1)

PR = 1;         % nearly-passive RIS  (no amplification budget)

% Transmit-power sweep [dBW]
p_min = -40; p_max = 40; p_step = 2;
[P_Watts, Ni] = available_power(p_min, p_max, p_step);

% Static power consumption of a NEARLY-PASSIVE RIS [dBm]  (Sec. V)
Po_dBm     = 40;   % all nodes except the RIS               (P0)
Pcn_dBm    = 0;    % per nearly-passive RIS element         (Pc,n^(p))
Po_ris_dBm = 20;   % other static sources at the RIS        (P0,RIS^(p))
Pc = static_power_consumption(Po_dBm, Pcn_dBm, Po_ris_dBm, N);

% Receiver noise   (PSD -174 dBm/Hz, noise figure 10 dB).
% A nearly-passive RIS adds no noise, so there is no sigma_RIS term.
noise_psd = -174; noise_figure = 10;
sigma_sq = noise_power(noise_psd, noise_figure, BW);

% Network geometry  (Sec. V)
r_cell    = 100;   % users within 100 m of the RIS
hbs       = 10;    % BS height  [m]
hris      = 15;    % RIS height [m]
hmax_ue   = 5;     % maximum user height [m]
ris_posix = 50;    % BS-to-RIS ground distance [m]

Ncarlo = 10;

%% ------------------------------------------------------------------ %%
%%  2. CHANNEL GENERATION                                              %%
%% ------------------------------------------------------------------ %%
G = zeros(NR, N, Ncarlo);
H = zeros(N,  K, Ncarlo);
for mc = 1:Ncarlo
    [G(:,:,mc), H(:,:,mc)] = generate_channels(r_cell, hbs, hris, hmax_ue, ris_posix, N, K, NR);
end
G = G(:,:,1);  H = H(:,:,1);   % use the first realization (see note in run_active_ris.m)

%% ------------------------------------------------------------------ %%
%%  3. PRE-ALLOCATION                                                  %%
%% ------------------------------------------------------------------ %%
SR1      = zeros(Ni,1);  GEE1      = zeros(Ni,1);
SR1_ropt = zeros(Ni,1);  GEE1_ropt = zeros(Ni,1);
SR1_eopt = zeros(Ni,1);  GEE1_eopt = zeros(Ni,1);
SR2_ropt = zeros(Ni,1);  GEE2_ropt = zeros(Ni,1);
SR2_eopt = zeros(Ni,1);  GEE2_eopt = zeros(Ni,1);
gamma_ropt = zeros(N,Ni); gamma_eopt = zeros(N,Ni);
X_ropt = zeros(N,N,Ni);   X_eopt = zeros(N,N,Ni);

opt_SR  = 0;   % maximize sum-rate
opt_GEE = 1;   % maximize GEE

% Feasible nearly-passive RIS start: unit global budget  ||gamma||^2 = N*PR.
phi   = exp(1i*2*pi*rand(N,1));
rho   = sqrt((N*PR) / norm(phi)^2);
gamma = rho .* phi;
X     = gamma * gamma';

%% ------------------------------------------------------------------ %%
%%  4. POWER SWEEP                                                     %%
%% ------------------------------------------------------------------ %%
for i = 1:Ni

    Pi_max = P_Watts(i);
    Pi     = (Pi_max / K) * ones(K,1);
    fprintf('\n=== Pmax = %g W (point %d/%d) ===\n', Pi_max, i, Ni);

    if i == 1
        P1_rin = Pi; P1_ein = Pi; P2_rin = Pi; P2_ein = Pi;
        gamma_rin = gamma; gamma_ein = gamma; X_rin = X; X_ein = X;
    end

    % (1) Uniform baseline
    C       = LMMSE_receiver_passive(Pi, G, H, gamma, sigma_sq, K, NR);
    sinr_p  = SINR_passive(Pi, C, G, H, gamma, sigma_sq, K);
    [~, SR1(i)] = data_rate(sinr_p, K, BW0);
    GEE1(i) = SR1(i) / (mu*sum(Pi) + Pc);

    % (2)-(3) Approach 1 - Algorithm 3
    [P1_rin, gamma_ropt(:,i), ~, SR1_ropt(i), GEE1_ropt(i)] = ...
        alt_opt1(P1_rin, Pi_max, G, H, gamma_rin, sigma_sq, K, NR, PR, BW0, mu, Pc, opt_SR);
    gamma_rin = gamma_ropt(:,i);

    [P1_ein, gamma_eopt(:,i), ~, SR1_eopt(i), GEE1_eopt(i)] = ...
        alt_opt1(P1_ein, Pi_max, G, H, gamma_ein, sigma_sq, K, NR, PR, BW0, mu, Pc, opt_GEE);
    gamma_ein = gamma_eopt(:,i);

    % (4)-(5) Approach 2 - Algorithm 6
    [P2_rin, X_ropt(:,:,i), SR2_ropt(i), GEE2_ropt(i)] = ...
        alt_opt2(P2_rin, Pi_max, G, H, X_rin, sigma_sq, K, NR, PR, BW0, mu, Pc, opt_SR);
    X_rin = X_ropt(:,:,i);

    [P2_ein, X_eopt(:,:,i), SR2_eopt(i), GEE2_eopt(i)] = ...
        alt_opt2(P2_ein, Pi_max, G, H, X_ein, sigma_sq, K, NR, PR, BW0, mu, Pc, opt_GEE);
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
title('Nearly-passive RIS: sum rate vs. transmit power');

figure;
plot(P_dBW, GEE1,'-ob', P_dBW, GEE1_ropt,'-dk', P_dBW, GEE2_ropt,'-og', ...
     P_dBW, GEE1_eopt,'-ro', P_dBW, GEE2_eopt,'-mo', 'LineWidth', 2);
grid on;
legend('Uniform (no opt.)', 'Approach 1 - SR opt', 'Approach 2 - SR opt', ...
       'Approach 1 - GEE opt', 'Approach 2 - GEE opt', 'Location', 'best');
xlabel('Maximum available power P_{max} [dBW]');
ylabel('Global energy efficiency [Mbit/J]');
title('Nearly-passive RIS: GEE vs. transmit power');
