function [btable,ttable] = load5MWDistributed

%% ----- Define Blade ----- %%
% Define constant parameters
Lb = 63;

% Load in distributed parameters
blade = readtable("Model Files/blade_distributed.txt");

% Define radial stations (minus hub)
blade.Rb = Lb * blade.BlFract;

% Define blade fundamental mode polynomials
F1_coeffs = [-2.2555, 4.7131, -3.2452, 1.7254, 0.0622, 0, 0];
E1_coeffs = [-0.6952, 2.376, -3.5772, 2.5337, 0.3627, 0, 0];

% Adjust coefficients to be function of x alone
idx = [6, 5, 4, 3, 2, 1, 0];
F1_coeffs = F1_coeffs./Lb.^idx;
E1_coeffs = E1_coeffs./Lb.^idx;

dF1_coeffs = polyder(F1_coeffs);
dE1_coeffs = polyder(E1_coeffs);

ddF1_coeffs = polyder(dF1_coeffs);
ddE1_coeffs = polyder(dE1_coeffs);

% Evaluate mode shapes and derivatives
blade.Fmode1 = polyval(F1_coeffs, blade.Rb);
blade.Emode1 = polyval(E1_coeffs, blade.Rb);

blade.dFmode1 = polyval(dF1_coeffs, blade.Rb);
blade.dEmode1 = polyval(dE1_coeffs, blade.Rb);

blade.ddFmode1 = polyval(ddF1_coeffs, blade.Rb);
blade.ddEmode1 = polyval(ddE1_coeffs, blade.Rb);

% % Compute in & out of plane mode shapes
% for i = 2:length(blade.Fmode1)
%     inner_flap(i) = trapz(blade.Rb(1:i), blade.ddFmode1(1:i).*cosd(blade.StrcTwist(1:i)));
%     inner_edge(i) = trapz(blade.Rb(1:i), blade.ddEmode1(1:i).*sind(blade.StrcTwist(1:i)));
% end
% 
% blade.

blade.EIo = blade.FlpStff .* cosd(blade.StrcTwst) + blade.EdgStff .* sind(blade.StrcTwst);


%% ----- Define Tower ----- %%
% Define constant parameters
Ht = 77.6;

% Read in distributed parameters
tower = readtable("Model Files/tower_distributed.txt");

% Define height stations
tower.Ht = tower.HtFract * Ht;

% Define fundamental mode
T1_coeffs = [0.4082, -1.5035, 1.8042, -0.8622, 1.1533, 0, 0];
T1_coeffs = T1_coeffs./Ht.^idx;

dT1_coeffs = polyder(T1_coeffs);

ddT1_coeffs = polyder(dT1_coeffs);

% Evaluate mode shape and derivatives
tower.mode1 = polyval(T1_coeffs, tower.Ht);
tower.dmode1 = polyval(dT1_coeffs, tower.Ht);
tower.ddmode1 = polyval(ddT1_coeffs, tower.Ht);

% figure
% plot(tower.Ht, tower.mode1)
% title('Fundamental Tower Mode')
% xlabel('Radius [m]')
% ylabel('Modal Displacement [m/m]')

btable = blade;
ttable = tower;

%% ----- Helper Functions ----- %%
% function val = innerIntegral(x, ddphi, theta)
%     idx = dsearchn(blade.Rb, x);
%     val = trapz(blade.Rb(1:idx), ddphi(1:idx).*cos(theta(1:idx)));
% end