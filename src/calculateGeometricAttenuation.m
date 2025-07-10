function attenuation = calculateGeometricAttenuation(d)
    % Constants
    c = 343; % Speed of sound in air at 20°C (in m/s)
    radius = d / 2; % Convert diameter to radius

    % Calculate geometric attenuation
    attenuation = 20*log10(1/(4*pi*radius));
end
