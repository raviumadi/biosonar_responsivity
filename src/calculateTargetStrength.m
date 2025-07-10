function TS = calculateTargetStrength(freqRange, diameter, distance, alpha)
% Calculates the acoustic Target Strength (TS) of a spherical reflector in air
% over a range of frequencies. This function is based on a simplified model
% and provides an approximation of TS.
%
% Inputs:
%   freqRange - A vector of frequencies (in Hz) at which TS will be calculated.
%               Example: linspace(100000, 10000, 100) for a downsweep from 100kHz to 10kHz.
%   diameter  - The diameter of the spherical reflector (in meters).
%   distance  - The distance from the reflector to the receiver (in meters).
%   alpha     - The absorption coefficient of air (in dB/m), which can depend
%               on frequency, humidity, and temperature.
%
% Output:
%   TS        - A vector of Target Strength values (in dB) corresponding to
%               each frequency in freqRange.
%
% Example usage:
%   freqRange = linspace(100000, 10000, 100); % Frequency range from 100kHz to 10kHz
%   diameter = 0.01; % Reflector diameter of 1 cm
%   distance = 1; % Distance of 1 meter
%   alpha = 0.1; % Absorption coefficient of 0.1 dB/m
%   TS = calculateTargetStrength(freqRange, diameter, distance, alpha);
%   plot(freqRange, TS);
%   xlabel('Frequency (Hz)');
%   ylabel('Target Strength (dB)');
%   title('Target Strength vs. Frequency');

% Constants
c = 343; % Speed of sound in air at 20°C (in m/s)
radius = diameter / 2; % Convert diameter to radius

% Calculate TS for each frequency
TS = 20*log10((2*pi*freqRange*radius)/c) - alpha*distance;

end
