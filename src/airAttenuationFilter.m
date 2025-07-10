function [ir, TL] = airAttenuationFilter(bl, fs, irpts, aatt)
% airAttenuationFilter Generates an impulse response with atmospheric and
% geometric attenuation.
%
% Inputs:
%   bl    - Calculated beam length minus physical and IO delays [meters]
%   fs    - Sampling frequency [Hz]
%   irpts - Length of the IR [samples]
%   aatt  - Additional attenuation due to object shape
%
% Outputs:
%   ir    - The impulse response with atmospheric and geometric attenuation
%   TL    - Geometric attenuation, aka TL.

% Define frequency bands and corresponding atmospheric attenuation values
f = [0 1e4 2e4 3e4 4e4 5e4 6e4 7e4 8e4 9e4 9.6e4] / (fs/2); % 10e4 11e4 12e4 13e4] / (fs/2);
dbmags = [0 -0.11 -0.4 -0.83 -1.31 -1.82 -2.33 -2.82 -3.3 -3.77 -4.06] *bl; % -4.2464 -4.7224 -5.2077 -5.7067] * bl;
mags = 10.^(dbmags ./ 20);

% Geometric attenuation calculation based on the beam length
scttr_coeff = 0.5;
b0 = 0.1; % Reference distance
p0 = 1; % Arbitrary sound pressure
% geom_atten = scttr_coeff * p0 * (b0 / bl)^2;


% Generate the impulse response with atmospheric attenuation
ir = fir2(irpts-1, f, mags);
attenuation = 2*calculateGeometricAttenuation(bl);
TL = max(ir) - 20e-6*10^(abs(attenuation)/20); % multiply this with the convolved signal and reduce by alpha

% Apply geometric and additional shape-based attenuation
% ir = ir .* geom_atten .* aatt;

end
