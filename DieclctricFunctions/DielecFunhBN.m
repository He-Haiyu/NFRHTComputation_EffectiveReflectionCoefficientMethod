function var = DielecFunhBN(omega,T)
% Ref: Near-Field heat transfer between graphene/hBN multilayers (Bo Zhao, 2017)
%% 
c = 3e10; % light, speed in cm/s

omg = omega/2/pi/c; % omg in cm^-1, omega in rad/s

%% perpendicular to optical axis (in plane)
epsiInf_perp = 4.87;
Gamma_perp = 5; % in cm^-1
OmegaTO_perp = 1370; % in cm^-1
OmegaLO_perp = 1610; % in cm^-1

epsi_perp = epsiInf_perp*(1 + (OmegaLO_perp^2 - OmegaTO_perp^2)./...
(OmegaTO_perp^2 - 1i*Gamma_perp*omg - omg.^2));

%% parallel to optical axis (out of plane)
epsiInf_para = 2.95;
Gamma_para = 4;
OmegaTO_para = 780;
OmegaLO_para = 830;

epsi_para = epsiInf_para*(1 + (OmegaLO_para^2 - OmegaTO_para^2)./...
(OmegaTO_para^2 - 1i*Gamma_para*omg - omg.^2));

%%
var = {epsi_perp,epsi_para};

end