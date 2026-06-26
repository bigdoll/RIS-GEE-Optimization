function [channel_G, channel_H] = generate_channels(r_cell, hbs, hris, hmax_ue, ris_posix, N, K, NR, nh, ng)
% GENERATE_CHANNELS  Geometry-based Rician channels for the RIS-aided uplink (G: RIS->BS, H: users->RIS).
%
%   function [channel_G, channel_H] = generate_channels(r_cell, hbs, hris, hmax_ue, ris_posix, N, K, NR, nh, ng)
%
% nh (optional): path-loss exponent of the user->RIS links  (default 4).
% ng (optional): path-loss exponent of the RIS->BS link     (default 4).
%   Distinct per-link exponents are used by the RIS-placement study (Fig. 7),
%   where the optimal RIS position depends on which link is more attenuated.
%   Calling with 8 arguments reproduces the single-exponent default (nh = ng = 4).
%
% Paper reference: System model, Sec. II / V  (per-link exponents: Fig. 7)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    if nargin < 9  || isempty(nh); nh = 4; end
    if nargin < 10 || isempty(ng); ng = 4; end

    % Cellular coverage r = 50m, D = 2r = 100m;
    r = r_cell; % radius of cell
    D =  2*r; % diameter of cell
    h_bs = hbs; % hieght of reciever (base_station)
    h_ris = hris; %1.5*h_bs; % height of RIS
    h_ue = hmax_ue*rand(K,1); % height of UES (hmax_ue = 10m)
    
    % position is defined as coordinates (horizontal_posix,vertical_posix)
    Rx = [0, h_bs]; % position of base station
    Ris = [ris_posix, h_ris]; % position of RIS
    d_ue = ris_posix + 10; % initial point for UEs (assume 10m from foot of RIS)
    Tx = [d_ue + (D-d_ue)*rand(K,1), h_ue]; % position of UEs

    % Path length from Transmitter (Tx) to (RIS) 
    dtx_ris = sqrt((Ris(:,1) - Tx(:,1)).^2 + (Ris(:,2) - Tx(:,2)).^2);

    % Path length from Receiver (Rx) to (RIS) 
    drx_ris = sqrt((Rx(:,1) - Ris(:,1)).^2 + (Rx(:,2) - Ris(:,2)).^2);

    % free space path-loss
    do = 35; % reference distance
    fo =  3.5e+9; % carrier frequency 3.5GHz
    c = 3e+8; % speed of light
    PLo = (4*pi*fo*do/c)^(-2); % free-space path loss at reference distance
    alpha_h = sqrt(2*PLo) ./ sqrt(1 + (dtx_ris/do).^(nh)); % path-loss coefficient between UE - RIS  (exponent nh)
    alpha_g = sqrt(2*PLo) ./ sqrt(1 + (drx_ris/do).^(ng)); % path-loss coefficient between RIS - Base Station (exponent ng)

    % Wireless channel for Tx - RIS ~ Rician fading
    K1 = 2;   K2 = 4;
    % H_LOS = exp(-1i*c_phase_h) .* alpha_h;
    H_NLOS = (randn(N, K) + 1i*randn(N, K)) / sqrt(2);
    G_NLOS = (randn(NR, N) + 1i*randn(NR, N)) / sqrt(2);
    channel_H=zeros(N,K);
    for k=1:K
           channel_H(:,k)=alpha_h(k)*(sqrt(K1/(K1+1))  + sqrt(1/(K1+1))*H_NLOS(:,k));
    end
    channel_G=alpha_g*(sqrt(K2/(K2+1))  + sqrt(1/(K2+1))*G_NLOS);
end
