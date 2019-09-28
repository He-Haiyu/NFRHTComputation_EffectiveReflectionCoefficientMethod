function var = DielecFunSi(omega,T)
% Ref: Handbook of Optical Constants of Solids, Edward D. Palik, ed. Academic Press, Boston, 1985 (from Github: nk)
% lambdan from 6e-4 microns to 2 microns; lambdak from 6e-4 to 1.17 microns
% omega from ~9e14 rad/s to ~3e18 rad/s
%%
c = 3e8; % light speed, in m/s
lambda = c ./ (omega/2/pi); % wavelength, in m
lambda = lambda * 1e6; % wavelength, in micron

%%
% n fitting
tem = xlsread('nk_Si_Palik_1985_0-2micron.xlsx',1);
lambdan = tem(2:end,1);
ns = tem(2:end,2);
[cfn,goodness] = fit(lambdan, ns, 'pchipinterp');
% k fitting
tem = xlsread('nk_Si_Palik_1985_0-2micron.xlsx',2);
lambdak = tem(2:end,1);
ks = tem(2:end,2);
[cfk,goodness] = fit(lambdak, ks, 'pchipinterp');

%%
% epsilon
epsi = (cfn(lambda) + cfk(lambda)*1i).^2;

var = {epsi,epsi};

end