function var = DielecFunSiO2(omega,T)
% Ref: Optical constants of silica glass from extreme ultraviolet to far infrared at near room temperature (Rei Kitamura, Laurent Pilon, 2007)
% for high frequencies, Im(epsilon) is neglegible. (wavelength < 7 μm, omega > 2.7e14 rad/s)
% fitted Sellmeier equation

%%
c = 3e8; % light speed, in m/s
lambda = c ./ (omega/2/pi); % wavelength, in m
lambda = lambda * 1e6; % wavelength, in μm

%% 
% 0.21 to 3.7 μm at 20 °C, not sure if can be used at high temperature
epsi_real = 1 + 0.6961663*lambda.^2./(lambda.^2-0.0684043^2) + ...
0.4079426*lambda.^2./(lambda.^2-0.1162414^2) + 0.8974794*lambda.^2./(lambda.^2-9.896161^2);
epsi_imag = 1e-10;

epsi = epsi_real + 1i*epsi_imag;

var = {epsi,epsi};

end
