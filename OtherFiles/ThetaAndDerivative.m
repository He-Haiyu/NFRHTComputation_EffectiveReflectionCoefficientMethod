clear;

% All temperatures interested in
Ts = [300,500,700,1000,1200];

%%
% Physical constants in SI unit
c0 = 299792458; % Speed of light in vacuum, m/s
ec = 1.60217646E-19; % Elementary charge, C
epsilon0 = 8.8541878176e-12; % Vacuum permitivity, F/m
hP = 6.626068e-34;  % Planck constant, J*s
hPb = 6.626068e-34/2/pi;  % Reduced Planck constant, J*s
kB = 1.3806503e-23; % Boltzmann constant, J/K
eV2cm = 8065.5443; % Conversion eV to 1/cm
cm2radps = 2*pi*c0*100; % Conversion from 1/cm to rad/s

%% frequency
Nfreq = 1000;
% omgeV = logspace(log10(0.0025),log10(0.32),Nfreq)'; % Energy, eV
% omg = omgeV*ec/hPb; % Angular frequency, rad/s
omg = logspace(log10(1e13),log10(5e14),Nfreq)'; % Angular frequency, rad/s
omgeV = omg*hPb/ec; % Energy, eV

%% Theta
Thetaomg = zeros(length(omg),length(Ts));

for i = 1:1:length(Ts)
    T = Ts(i);
    Thetaomg(:,i) = hPb*omg ./ (exp(hPb*omg/kB/T) - 1);
end

%% dTheta/dT

Thetaderivative = zeros(length(omg),length(Ts));
dT = 0.0001;
for j = 1:1:length(Ts)
    T = Ts(j);
    Th = T + dT;
    Tl = T - dT;
    Thetaderivative(:,j) = (hPb*omg ./ (exp(hPb*omg/kB/Th) - 1) - hPb*omg ./ (exp(hPb*omg/kB/Tl) - 1))/(2*dT);
end

%% Initialize figures
% Theta
fignm = ['Theta:' ' T=' num2str(T)];
h1 = figure('Name',fignm); ah1 = gca;
set(ah1,'Fontname','times','Fontsize',12,'Yscale','log','Box','on');
xlabel('\omega (rad/s)');
ylabel('\Theta_\omega (W m^-^2 rad^-^1 s)');

axis([omg(1),omg(end),-inf,inf]);
hold all;

clear fignm figtitle;

% dTheta/dT
fignm = ['dTheta/dT: ' 'T=' num2str(T) ' dT=' num2str(dT)];
h2 = figure('Name',fignm); ah2 = gca;
set(ah2,'Fontname','times','Fontsize',12,'Yscale','log','Box','on');
xlabel('\omega (rad/s)');
ylabel('d\Theta_\omega / dT (W m^-^2 rad^-^1 s K^-^1)');

axis([omg(1),omg(end),-inf,inf]);
hold all;

clear fignm figtitle;

%% Plot
for i = 1:length(Ts)
    plot(ah1,omg,Thetaomg(:,i),'LineWidth',1.5);
    plot(ah2,omg,Thetaderivative(:,i),'LineWidth',1.5);
end



