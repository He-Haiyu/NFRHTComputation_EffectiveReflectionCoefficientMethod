%% Constants
% Physical constants in SI unit
c0 = 299792458; % Speed of light in vacuum, m/s
ec = 1.60217646E-19; % Elementary charge, C
epsilon0 = 8.8541878176e-12; % Vacuum permitivity, F/m
hP = 6.626068e-34;  % Planck constant, J*s
hPb = 6.626068e-34/2/pi;  % Reduced Planck constant, J*s
kB = 1.3806503e-23; % Boltzmann constant, J/K
eV2cm = 8065.5443; % Conversion eV to 1/cm
cm2radps = 2*pi*c0*100; % Conversion from 1/cm to rad/s

%% Initiation
LayerNames = {'VO2','VO2'}; % from top to bottom
LayerThicknesses = [1,1];
T = 330;

N = length(LayerNames);
d = LayerThicknesses; % for simplicity
%% Omega
% 
Nfreq = 1000;
omg = logspace(log10(1e13),log10(3e14),Nfreq)';

%% Dielectric Properties
epslver = zeros(length(omg),N); % perpendicular to optical axis (in plane)
epslpara = zeros(length(omg),N); % parallel to optical axis (out of plane)

for i = 1:N
    epslif = str2func(['DielecFun',LayerNames{i}]);
    epsli = epslif(omg,T);
    epslver(:,i) = epsli{1};
    epslpara(:,i) = epsli{2};
end   

%% 
% Mean energy of harmonic oscillator without zero-point hPb*omg/2
Thetaomg = hPb*omg ./ (exp(hPb*omg/kB/T) - 1);

%% MAIN PART
% Numerical error control
rTol = 1e-7;
aTol = 1e-7;
% Initialize q arry
qs = zeros(length(omg));
qp = zeros(length(omg));
q = qs + qp;


for j = 1:length(omg) % one j for one omega value
    
    omgj = omg(j);
    k0j = omgj/c0;
    
    epslverj = epslver(j,:);
    epslparaj = epslpara(j,:);

    Thetaomgj = Thetaomg(j);
        
    % Wavevector component perpendicular to interface
    gm0 = @(x) sqrt(k0j^2-x.^2); 
    
    gms = cell(1,N);
    gmp = cell(1,N);
    for i=1:1:N  
        gms{i}= @(x) checksqrt(epslverj(i)*k0j^2-x.^2);
        gmp{i}= @(x) checksqrt(epslverj(i)*k0j^2-epslverj(i)*x.^2/epslparaj(i)); % 1/m
    end

    
    % Fresnel reflection coefficient: interface
    % ri stands for the ith interface, not including the media-vacuum interface
    r01s = @(x) (gm0(x)-gms{1}(x)) ./ (gm0(x)+gms{1}(x)); % s-polarization
    r01p = @(x) (epslverj(1).*gm0(x)-gmp{1}(x)) ./ (epslverj(1).*gm0(x)+gmp{1}(x)); % p-polarization

    
    rs = cell(1,N-1);
    rp = cell(1,N-1);

    
    for i=1:1:N-1
        rs{i} = @(x) (gms{i}(x)-gms{i+1}(x)) ./ (gms{i}(x)+gms{i+1}(x));
        rp{i} = @(x) (epslverj(i+1).*gmp{i}(x)-epslverj(i).*gmp{i+1}(x)) ./ (epslverj(i+1).*gmp{i}(x)+epslverj(i).*gmp{i+1}(x));
    end
    
    % Fresnel reflection coefficient: multilayer system
    % R for generalized fresnel reflection coefficient
    
    Rs = cell(1,N-1);
    Rp = cell(1,N-1);
   
    Rs{N-1} = rs{N-1};
    Rp{N-1} = rp{N-1};
    
    for i=1:1:N-2
        %interface between layer N-1-i and N-i
        Rs{N-1-i} = @(x) (rs{N-1-i}(x)+Rs{N-i}(x).*exp(2j*gms{N-i}(x)*d(N-i)))...
            ./ (1+rs{N-1-i}(x).*Rs{N-i}(x).*exp(2j*gms{N-i}(x)*d(N-i))); 
        Rp{N-1-i} = @(x) (rp{N-1-i}(x)+Rp{N-i}(x).*exp(2j*gmp{N-i}(x)*d(N-i)))...
            ./ (1+rp{N-1-i}(x).*Rp{N-i}(x).*exp(2j*gmp{N-i}(x)*d(N-i)));
    end
    
    Rs0 = @(x) (r01s(x)+Rs{1}(x).*exp(2j*gms{1}(x)*d(1))) ./ (1+r01s(x).*Rs{1}(x).*exp(2j*gms{1}(x)*d(1)));
    Rp0 = @(x) (r01p(x)+Rp{1}(x).*exp(2j*gmp{1}(x)*d(1))) ./ (1+r01p(x).*Rp{1}(x).*exp(2j*gmp{1}(x)*d(1)));
    
    % Integral Term
    ITs = @(x) (1 - abs(Rs0(x)).^2);
    ITp = @(x) (1 - abs(Rp0(x)).^2);
    
    % Far field energy flux
    qs(j) = integral(@(x) ITs(x).*x,0,k0j,'RelTol',rTol,'AbsTol',aTol)*Thetaomgj/2/pi^2;
    qp(j) = integral(@(x) ITp(x).*x,0,k0j,'RelTol',rTol,'AbsTol',aTol)*Thetaomgj/2/pi^2;
    
    % inflector
    if mod(j,100) == 0
        j
    end
end

% q
q = qs + qp;


%% Plot
% Initialize figure
fignm = ['SpecQ'];
h = figure('Name',fignm); ah = gca;
set(ah,'Fontname','times','Fontsize',12,'Yscale','log');
xlabel('\omega (rad/s)');
ylabel('Q_\omega (W m^-^2 rad^-^1 s)');
axis([omg(1),omg(end),-inf,inf]);
hold all;

clear fignm figtitle;
% plot
plot(ah,omg,abs(q),'LineWidth',1.5);

