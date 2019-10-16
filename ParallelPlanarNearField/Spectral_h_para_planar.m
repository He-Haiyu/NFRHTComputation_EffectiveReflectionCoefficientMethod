%%
% Clear the enviroment
% close all;
clear;clc;format long

%% INI, at least two layers for each side
MaterialB = {'BN','BN'}; % from the gap size
MaterialT = {'VO2','VO2'}; % from the gap size
ThicknessB = [5e-7,1];
ThicknessT = [0.5e-6,1];

gp = 1e-7; % gap size

% Temperature
T = 330;
% T difference
dT = 0.001;

% Suppose top structure has higher temperature
Tt = T + dT; % top structure temperature
Tb = T - dT; % bottom structure temperature

Nb = length(MaterialB);
Nt = length(MaterialT);

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

%% Omega and k
Nfreq = 1000;
% omgeV = logspace(log10(0.0025),log10(0.32),Nfreq)'; % Energy, eV
% omg = omgeV*ec/hPb; % Angular frequency, rad/s

omg = logspace(log10(1e13),log10(3e14),Nfreq)'; % Angular frequency, rad/s
omgeV = omg*hPb/ec; % Energy, eV

% nv = omg/2/pi; % Frequency, Hz
% lmd = c0./nv; % Wavelength, m
% k = 1./lmd/100; % Wavenumber, 1/cm

%% Physical properties
% Bottom
epslbver = zeros(length(omg),Nb); % Ordinary dielectric function, perpendicular to the optical axis, parallel to the surface plane
epslbpara = zeros(length(omg),Nb); % Extraordinary dielectric function, parallel to the optical axis, perpendicular to the surface plane

for i = 1:Nb
    epslbif = str2func(['DielecFun',MaterialB{i}]);
    epslbi = epslbif(omg,Tb);
    epslbver(:,i) = epslbi{1};
    epslbpara(:,i) = epslbi{2};
end    

% Top
epsltver = zeros(length(omg),Nt); % Ordinary dielectric function, perpendicular to the optical axis, parallel to the surface plane
epsltpara = zeros(length(omg),Nt); % Extraordinary dielectric function, parallel to the optical axis, perpendicular to the surface plane

for i = 1:Nt
    epsltif = str2func(['DielecFun',MaterialT{i}]);
    epslti = epsltif(omg,Tt);
    epsltver(:,i) = epslti{1};
    epsltpara(:,i) = epslti{2};
end 

db = ThicknessB;
dt = ThicknessT;

%%
% Mean energy of harmonic oscillator without zero-point hPb*omg/2
ThetaomgTb = hPb*omg ./ (exp(hPb*omg/kB/Tb) - 1);
ThetaomgTt = hPb*omg ./ (exp(hPb*omg/kB/Tt) - 1);
dThetaomgbt = ThetaomgTt - ThetaomgTb;

%% 
% Numerical error control
rTol = 1e-7;
aTol = 1e-7;

for j = 1:length(omg)
    omgj = omg(j);
    k0j = omgj/c0;
    
    epslbverj = epslbver(j,:);
    epslbparaj = epslbpara(j,:);
    epsltverj = epsltver(j,:);
    epsltparaj = epsltpara(j,:);
     
    dThetaomgbtj = dThetaomgbt(j);
        
    % Wavevector component perpendicular to interface
    gm0 = @(x) sqrt(k0j^2-x.^2); % 1/m
    gmbs = cell(1,Nb);
    gmbp = cell(1,Nb);
    gmts = cell(1,Nt);
    gmtp = cell(1,Nt);
    for b=1:1:Nb
        gmbs{b}= @(x) checksqrt(epslbverj(b)*k0j^2-x.^2);
        gmbp{b}= @(x) checksqrt(epslbverj(b)*k0j^2-epslbverj(b)*x.^2/epslbparaj(b)); % 1/m
    end
    
    for t=1:1:Nt
        gmts{t}= @(x) checksqrt(epsltverj(t)*k0j^2-x.^2);
        gmtp{t}= @(x) checksqrt(epsltverj(t)*k0j^2-epsltverj(t)*x.^2/epsltparaj(t)); % 1/m     
    end
    
    % Fresnel reflection coefficient: interface
    % ri stands for the ith interface, not including the media-vacuum interface
    r0b1s = @(x) (gm0(x)-gmbs{1}(x)) ./ (gm0(x)+gmbs{1}(x)); % s-polarization
    r0b1p = @(x) (epslbverj(1).*gm0(x)-gmbp{1}(x)) ./ (epslbverj(1).*gm0(x)+gmbp{1}(x)); % p-polarization
    r0t1s = @(x) (gm0(x)-gmts{1}(x)) ./ (gm0(x)+gmts{1}(x)); % s-polarization
    r0t1p = @(x) (epsltverj(1).*gm0(x)-gmtp{1}(x)) ./ (epsltverj(1).*gm0(x)+gmtp{1}(x)); % p-polarization
    
    rbs = cell(1,Nb-1);
    rbp = cell(1,Nb-1);
    rts = cell(1,Nt-1);
    rtp = cell(1,Nt-1);
    
    for i=1:1:Nb-1
        rbs{i} = @(x) (gmbs{i}(x)-gmbs{i+1}(x)) ./ (gmbs{i}(x)+gmbs{i+1}(x));
        rbp{i} = @(x) (epslbverj(i+1).*gmbp{i}(x)-epslbverj(i).*gmbp{i+1}(x)) ./ (epslbverj(i+1).*gmbp{i}(x)+epslbverj(i).*gmbp{i+1}(x));
    end
    
    for i=1:1:Nt-1
        rts{i} = @(x) (gmts{i}(x)-gmts{i+1}(x)) ./ (gmts{i}(x)+gmts{i+1}(x));
        rtp{i} = @(x) (epsltverj(i+1).*gmtp{i}(x)-epsltverj(i).*gmtp{i+1}(x)) ./ (epsltverj(i+1).*gmtp{i}(x)+epsltverj(i).*gmtp{i+1}(x));
    end
    
    % Fresnel reflection coefficient: multilayer system
    % R for generalized fresnel reflection coefficient
    
    Rbs = cell(1,Nb-1);
    Rbp = cell(1,Nb-1);
    Rts = cell(1,Nt-1);
    Rtp = cell(1,Nt-1);    
    
    
    % Bottom
    Rbs{Nb-1} = rbs{Nb-1};
    Rbp{Nb-1} = rbp{Nb-1};
    
    for i=1:1:Nb-2
        %interface between layer Nb-1-i and Nb-i
        Rbs{Nb-1-i} = @(x) (rbs{Nb-1-i}(x)+Rbs{Nb-i}(x).*exp(2j*gmbs{Nb-i}(x)*db(Nb-i)))...
            ./ (1+rbs{Nb-1-i}(x).*Rbs{Nb-i}(x).*exp(2j*gmbs{Nb-i}(x)*db(Nb-i))); 
        Rbp{Nb-1-i} = @(x) (rbp{Nb-1-i}(x)+Rbp{Nb-i}(x).*exp(2j*gmbp{Nb-i}(x)*db(Nb-i)))...
            ./ (1+rbp{Nb-1-i}(x).*Rbp{Nb-i}(x).*exp(2j*gmbp{Nb-i}(x)*db(Nb-i)));
    end
    
    Rbs0 = @(x) (r0b1s(x)+Rbs{1}(x).*exp(2j*gmbs{1}(x)*db(1))) ./ (1+r0b1s(x).*Rbs{1}(x).*exp(2j*gmbs{1}(x)*db(1)));
    Rbp0 = @(x) (r0b1p(x)+Rbp{1}(x).*exp(2j*gmbp{1}(x)*db(1))) ./ (1+r0b1p(x).*Rbp{1}(x).*exp(2j*gmbp{1}(x)*db(1)));    
    
    
    %Top system
    Rts{Nt-1} = rts{Nt-1};
    Rtp{Nt-1} = rtp{Nt-1};
    
    for i=1:1:Nt-2
        %interface between layer Nt-1-i and Nt-i
        Rts{Nt-1-i} = @(x) (rts{Nt-1-i}(x)+Rts{Nt-i}(x).*exp(2j*gmts{Nt-i}(x)*dt(Nt-i)))...
            ./ (1+rts{Nt-1-i}(x).*Rts{Nt-i}(x).*exp(2j*gmts{Nt-i}(x)*dt(Nt-i))); 
        Rtp{Nt-1-i} = @(x) (rtp{Nt-1-i}(x)+Rtp{Nt-i}(x).*exp(2j*gmtp{Nt-i}(x)*dt(Nt-i)))...
            ./ (1+rtp{Nt-1-i}(x).*Rtp{Nt-i}(x).*exp(2j*gmtp{Nt-i}(x)*dt(Nt-i)));
    end
    
    Rts0 = @(x) (r0t1s(x)+Rts{1}(x).*exp(2j*gmts{1}(x)*dt(1))) ./ (1+r0t1s(x).*Rts{1}(x).*exp(2j*gmts{1}(x)*dt(1)));
    Rtp0 = @(x) (r0t1p(x)+Rtp{1}(x).*exp(2j*gmtp{1}(x)*dt(1))) ./ (1+r0t1p(x).*Rtp{1}(x).*exp(2j*gmtp{1}(x)*dt(1)));

	% Transmission coefficient
    % Propagating
    Tprs = @(x) (1-abs(Rbs0(x)).^2).*(1-abs(Rts0(x)).^2) ./ abs(1-Rbs0(x).*Rts0(x).*exp(2j*gm0(x)*gp)).^2;
    Tprp = @(x) (1-abs(Rbp0(x)).^2).*(1-abs(Rtp0(x)).^2) ./ abs(1-Rbp0(x).*Rtp0(x).*exp(2j*gm0(x)*gp)).^2;
    % Evanescent
    Tevs = @(x) 4*imag(Rbs0(x)).*imag(Rts0(x)).*exp(-2*abs(gm0(x))*gp) ./ abs(1-Rbs0(x).*Rts0(x).*exp(2j*gm0(x)*gp)).^2;
    Tevp = @(x) 4*imag(Rbp0(x)).*imag(Rtp0(x)).*exp(-2*abs(gm0(x))*gp) ./ abs(1-Rbp0(x).*Rtp0(x).*exp(2j*gm0(x)*gp)).^2;	
    
    %Energy flux due to omgj
    qpromgs(j) = integral(@(x) Tprs(x).*x,0,k0j,'RelTol',rTol,'AbsTol',aTol)*dThetaomgbtj/4/pi^2;
    qpromgp(j) = integral(@(x) Tprp(x).*x,0,k0j,'RelTol',rTol,'AbsTol',aTol)*dThetaomgbtj/4/pi^2;
    qevomgs(j) = integral(@(x) Tevs(x).*x,k0j,k0j+8*pi/gp,'RelTol',rTol,'AbsTol',aTol)*dThetaomgbtj/4/pi^2;
    qevomgp(j) = integral(@(x) Tevp(x).*x,k0j,k0j+8*pi/gp,'RelTol',rTol,'AbsTol',aTol)*dThetaomgbtj/4/pi^2;
    
    if mod(j,100) == 0 % flag
        j
    end
end

% Conductance
hpromgp = qpromgp/(2*dT);
hpromgs = qpromgs/(2*dT);
hevomgp = qevomgp/(2*dT);
hevomgs = qevomgs/(2*dT);

h = hpromgp+hpromgs+hevomgp+hevomgs;

hprs = trapz(omg,hpromgs);
hprp = trapz(omg,hpromgp);
hevs = trapz(omg,hevomgs);
hevp = trapz(omg,hevomgp);

% print heat fluxes
hunit = ' W/(m^2 rad/s K)';
disp(['In: ',hunit]);
disp(['h = ',num2str(hprs+hprp+hevs+hevp),'hs = ',num2str(hprs+hevs),', hp = ',num2str(hprp+hevp)]);
disp(['hpr = ',num2str(hprs+hprp),', hprs = ',num2str(hprs),', hprp = ',num2str(hprp)]);
disp(['hev = ',num2str(hevs+hevp),', hevs = ',num2str(hevs),', hevp = ',num2str(hevp)]);

%% Initialize figure
%Spectral_h
fignm = ['Spech: ' 'gp=' num2str(gp) ' T=' num2str(T) ', dT=' num2str(dT)];
h3 = figure('Name',fignm); ah = gca;
set(ah,'Fontname','times','Fontsize',12,'Yscale','log','Box','on');
xlabel('\omega (rad/s)');
ylabel('h_\omega (W m^-^2 rad^-^1 s K^-^1)');

axis([omg(1),omg(end),-inf,inf]);
hold all;

clear fignm figtitle;


%% Plot
plot(ah,omg,abs(h),'LineWidth',1.5);
plot(ah,omg,abs(hpromgs+hpromgp),'LineWidth',1.5);
plot(ah,omg,abs(hevomgs+hevomgp),'LineWidth',1.5);
legend('Total','Propagating','Evanescent');
