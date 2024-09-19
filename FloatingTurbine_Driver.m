clear all; close all; clc;

global btable ttable

%% ----- Load in Tower & Blade Information ----- %%
[btable, ttable] = load5MWDistributed;
blade = readmatrix('OpenFAST\5MW_Baseline\AeroData\NREL_5MW_Blade.csv');

%% ----- Load in OpenFAST Results ----- %%
% Steady 13 m/s simulation
% load('C:\Umaine Google Sync\Masters Working Folder\1 - OpenFAST\Simulations\5MW_OC4Semi_TrinityDOF\SteadyWind_13mps_NoWave\5MW_OC4Semi_TrinityDOF_FAST_Results.mat')

% Free decay (no aerodynamics)
% load('C:\Umaine Google Sync\Masters Working Folder\1 - OpenFAST\Simulations\5MW_OC4Semi_TrinityDOF\FD_Pitch_NoAero\5MW_OC4Semi_TrinityDOF_FAST_Results.mat');

% Turbulent 13 m/s simulation
load('OpenFAST\Simulations\5MW_OC4Semi_WSt_WavesWN\Turbulent_13mps_Monhegan_NoWave\5MW_OC4Semi_WSt_WavesWN_FAST_Results.mat');

rotAz = sim_results.Azimuth;

%% ----- Simulation Setup ----- %%
% Initial conditions
x0 = zeros(12,1);
x = x0;
% x(1) = 0.005;

% Time vector
dt = 0.01;
t = [0:dt:60];
N = length(t);

% Storage Variables
xstore = zeros(length(x0), length(t));

% Input structure
F = struct();

%% ----- Run Standard Simulation ----- %%
flaplog = zeros(3,length(t));
xdot = zeros(6,1);

Q = 10^7*eye(length(x));
R = 10^-6*eye(2);

tic
for i = 1:length(t)

    % Update azimuth measurement
    F.azimuth = deg2rad(rotAz(i));
    F.omega = sim_results.RotSpeed(i) * (pi/30);
    RPM = sim_results.RotSpeed(i);
    HH_wind = sim_results.Wind1VelX(i);

    % Project wind onto blades
    azimuth = deg2rad(sim_results.Azimuth(i));

    zblade_1 = 90 + cos(azimuth)*blade(:,1);
    zblade_2 = 90 + cos(azimuth + deg2rad(120))*blade(:,1);
    zblade_3 = 90 + cos(azimuth + deg2rad(240))*blade(:,1);

    ublade_1 = powerLaw(HH_wind,90,zblade_1,0.2);
    ublade_2 = powerLaw(HH_wind,90,zblade_2,0.2);
    ublade_3 = powerLaw(HH_wind,90,zblade_3,0.2);

    wind_vector = [ublade_1,ublade_2,ublade_3];

    % Get blade pitch commands
    blade_pitch = (pi/180)*[sim_results.BldPitch1(i), sim_results.BldPitch2(i), sim_results.BldPitch3(i)];

    % Get loads
    aero_forces = RT_BEM_Ning_v2(wind_vector,RPM,blade_pitch);
    % hydro_forces = RT_Hydro(dt, x, xdot);
    hydro_forces.PtfmPitch = 0;

    % Applied forces
    F.forces = zeros(6,1);
    F.forces(1) = sum(aero_forces.Thrust)*87.6 + hydro_forces.PtfmPitch;
    F.forces(2) = sum(aero_forces.Thrust);
    F.forces(3:5) = aero_forces.FlapwiseForce;

    % Update state vector
    xp = RK4(dt, x, @FloatingTrinityRHS, F);

    % Store state vector
    xdot = (xp - x)/dt;
    xstore(:,i) = xp;
    x = xp;
end
tottime = toc;
fprintf('CPU Time Ratio: %0.3g \n', tottime/(max(t)-min(t)));

%% ----- Plot Results ----- %%
figure
gca; hold on; box on;
plot(t, xstore(1,:),'DisplayName', 'Linear')
plot(sim_results.Time,sim_results.PtfmPitch*(pi/180), 'DisplayName', 'FAST');
xlabel('Time [s]')
xlim([min(t),max(t)])
ylabel('Displacement [rads]')
title('Platform Pitch Displacement')
legend

figure
gca; hold on; box on;
plot(t, xstore(2,:),'DisplayName', 'Linear')
plot(sim_results.Time,sim_results.TTDspFA, 'DisplayName', 'FAST');
xlabel('Time [s]')
xlim([min(t),max(t)])
ylabel('Displacement [m]')
title('Tower Top Displacement')
legend

figure
gca; hold on; box on;
plot(t, xstore(3,:),'DisplayName', 'Linear')
plot(sim_results.Time,sim_results.Q_B1F1, 'DisplayName', 'FAST');
xlabel('Time [s]')
xlim([min(t),max(t)])
ylabel('Displacement [m]')
title('Blade 1 Flapwise Modal Displacement')
legend

%% ---------- HELPER FUNCTIONS ---------- %%
function u = powerLaw(uref,zref,z,alpha)
    u = uref.*(z./zref).^alpha;
end