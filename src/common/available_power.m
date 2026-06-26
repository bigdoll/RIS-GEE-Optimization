function [P_Watts, Ni] = available_power(p_min, p_max, p_step)
% AVAILABLE_POWER  Build the transmit-power sweep: convert a dB gain range to linear Watts.
%
%   function [P_Watts, Ni] = available_power(p_min, p_max, p_step)
%
% Paper reference: System model, Sec. V
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.
    % Available power 
    P_dB = p_min:p_step:p_max; % Power gain range
    Ni = length(P_dB); % length of power array
    P_Watts = 10.^(P_dB/10); % power conversion from dB (dBW or dBm) to watts         
end
