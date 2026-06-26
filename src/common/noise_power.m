function sigma_sq = noise_power(noise_psd, noise_figure, BW)
% NOISE_POWER  Thermal-noise power (variance) from PSD, noise figure and bandwidth.
%
%   function sigma_sq = noise_power(noise_psd, noise_figure, BW)
%
% Paper reference: System model, Sec. II
% Part of: RIS-GEE-Optimization  (github repo). See docs/equation_mapping.md.

    % Noise Figure
    Fn_dB = noise_figure;
    Fn = 10^(0.1*Fn_dB);
    
    % Noise Power Spectal Density
    No_dBm = noise_psd;
    No = 10^((No_dBm-30)/10);
    
    % Noise Power or Variance
    sigma_sq = Fn * No * BW;
    
end
