function var = DielecFunGr(omega,T)
% Ref: Graphene-assisted near-field radiative thermal rectifier based on
% phase transition of vanadium dioxide (VO2) (Zheng Z., 2017)

mueV = 0; %eV %Chemical Potential

%Constants
hbar = 6.626068e-34/2/pi;%J*s %Reduced Planck constant
tau = 1E-13;%s %relaxation
tg = 0.34E-9;%m %graphene sheet thickness
kb = 1.3806503e-23;%J/K %Boltzmann Constant
e = 1.60217646E-19;%C %elementary charge
epslv= 8.854187817E-12;%F/m %vacuum permittivity

mu = mueV*e; %J %Chemical Potential

% G = @(x) sinh(x/kb/T)./(cosh(x/kb/T)+cosh(mu/kb/T));
G = @(x) (1-exp(-2*x))./(1+exp(-2*x)+2*cosh(mu/kb/T)*exp(-x));
% G = @(x) (exp(x)-exp(-x))./(exp(x)+exp(-x)+2*cosh(mu/kb/T));

SigmaDrudefun = @(ome) 1i./(ome+1i/tau) * e^2/pi/hbar^2 * 2*kb*T*log(2*cosh(mu/2/kb/T));

SigmaInter = zeros(length(omega),1);
% etaInt = zeros(1:length(omega));
for j = 1:length(omega)
    omegaj = omega(j);
    etafunj = @(x) (G(x)-G(hbar*omegaj/2)) ./ ((hbar*omegaj)^2 - 4*x.^2);
    etaIntj1 = integral(etafunj,0,1e-14);
    etaIntj2 = integral(etafunj,1e-14,1e-12);
    etaIntj3 = integral(etafunj,1e-12,inf);
    etaIntj = etaIntj1+etaIntj2+etaIntj3;
    SigmaInter(j) = e^2/4/hbar* (G(hbar*omegaj/2)+1i*4*hbar*omegaj/pi*etaIntj);
end

SigmaDrude = SigmaDrudefun(omega);

Sigma = SigmaInter + SigmaDrude;

epsl = 1+1i*Sigma./(omega*tg*epslv);

var = {epsl,epsl};

end