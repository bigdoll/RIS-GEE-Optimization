function [Pc] = static_power_consumption(Po_dBm, Pcn_dBm, Po_ris_dBm, N)
% STATIC_POWER_CONSUMPTION  Static circuit power  Pc = P0 + N*Pc,n + P0,RIS.
%
%   function [Pc] = static_power_consumption(Po_dBm, Pcn_dBm, Po_ris_dBm, N)
%
% Paper reference: Eq. (10)
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    % Static power consumption
    % Po_dBm total static power consumption of all other nodes of the system except the RIS in dBm
    % Pcn_dBm static power consumed by each RIS element in dBm
    % Po_ris_dBm all other static power consumption source of the RIS in dBm
    % N:  number of RIS reflecting elements
    Po = 10^((Po_dBm-30)/10); % conversion to Watts
    Pcn = 10^((Pcn_dBm-30)/10); % conversion to Watts
    Po_ris = 10^((Po_ris_dBm-30)/10); % conversion to Watts
    Pc = Po + N*Pcn + Po_ris; % total power consumed in the system
end
