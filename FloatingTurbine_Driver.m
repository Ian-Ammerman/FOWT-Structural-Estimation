clear all; close all; clc;

global btable ttable

%% ----- Load in Tower & Blade Information ----- %%
[btable, ttable] = load5MWOC4SemiDistributed;
blade = readmatrix('OpenFAST\5MW_Baseline\AeroData\NREL_5MW_Blade.csv');

%% ----- Load in OpenFAST Results ----- %%
% Turbulent 13 m/s simulation (FIXED-BOTTOM)
% load('OpenFAST\Simulations\5MW_Baseline_Fixed\Turbulent_13mps_Monhegan\5MW_Baseline_Fixed_FAST_Results.mat');

% Steady 13 m/s simulation
% load('C:\Umaine Google Sync\Masters Working Folder\1 - OpenFAST\Simulations\5MW_OC4Semi_TrinityDOF\SteadyWind_13mps_NoWave\5MW_OC4Semi_TrinityDOF_FAST_Results.mat')

% Free decay (no aerodynamics)
% load('C:\Umaine Google Sync\Masters Working Folder\1 - OpenFAST\Simulations\5MW_OC4Semi_TrinityDOF\FD_Pitch_NoAero\5MW_OC4Semi_TrinityDOF_FAST_Results.mat');

% Turbulent 13 m/s simulation
load('OpenFAST\Simulations\5MW_OC4Semi_WSt_WavesWN\Turbulent_13mps_Monhegan_NoWave\5MW_OC4Semi_WSt_WavesWN_FAST_Results.mat');

% Turbulent 13 m/s simulation w/ Waves
% load('OpenFAST\Simulations\5MW_OC4Semi_WSt_WavesWN\Turbulent_13mps_Monhegan_Hs3d2_Tp12\5MW_OC4Semi_WSt_WavesWN_FAST_Results.mat');

rotAz = sim_results.Azimuth;

%% ----- Simulation Setup ----- %%
% Initial conditions
x0 = zeros(12,1);
x = x0;
% x(1) = 0.005;

% Time vector
dt = 0.01;
t = [0:dt:250];
N = length(t);

% Storage Variables
xstore = zeros(length(x0), length(t));

% Input structure
F = struct();

%% ----- Run Standard Simulation ----- %%
flaplog = zeros(3,length(t));
xdot = zeros(12,1);

Q = 10^7*eye(length(x));
R = 10^-6*eye(2);

tic
for i = 1:length(t)

    % Update azimuth measurement
    F.azimuth = deg2rad(rotAz(i));
    F.omega = sim_results.RotSpeed(i) * (pi/30);
    RPM = sim_results.RotSpeed(i);
    dqb = [x(9),x(10),x(11)];
    Uinf = sim_results.Wind1VelX(i);
    HH_wind = Uinf - x(8) - 87.6*x(7);
    kf_measurements = [sim_results.QD_P(i);
                       sim_results.QD_TFA1(i)];

    % Project wind onto blades
    wind_vector = projectWind(HH_wind, F.azimuth, blade, 90, 0.18);

    % Get blade pitch commands
    blade_pitch = (pi/180)*[sim_results.BldPitch1(i), sim_results.BldPitch2(i), sim_results.BldPitch3(i)];

    % Get loads
    aero_forces = RT_BEM_Ning_v2(wind_vector,dqb,RPM,blade_pitch);
    Ftwr_aero = TwrAeroRT(Uinf,x(7),0.18);
    hydro_forces = RT_Hydro(dt, x, xdot);
    % hydro_forces.PtfmPitch = -sim_results.HydroMyi(i);
    % hydro_forces.PtfmPitch = 0;

    % Applied forces
    F.forces = zeros(6,1);
    F.forces(1) = sum(aero_forces.Thrust)*87.6 + hydro_forces.PtfmPitch;
    F.forces(2) = sum(aero_forces.Thrust) + Ftwr_aero;
    F.forces(3:5) = aero_forces.FlapwiseForce;
    F.forces(6) = 0;
    % F.forces(3:5)

    % Update state vector
    xp = RK4(dt, x, @FloatingTrinityRHS, F);
    % xp = SRCUKF(dt, x, @FloatingTrinityRHS, @FloatingTrinityMeasurement, Q, R, kf_measurements, F);

    % Store state vector
    xdot = (xp - x)/dt;
    xstore(:,i) = xp;
    x = xp;

    % Print elapsed time
    if rem(t(i),10) == 0
        fprintf('Elapsed Time: %0.3g seconds \n',t(i))
    end
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

function wind_vector = projectWind(HHwind, RotAz, Blade, HubHt, alpha)
    zblade_1 = 90 + cos(RotAz)*Blade(:,1);
    zblade_2 = 90 + cos(RotAz + deg2rad(120))*Blade(:,1);
    zblade_3 = 90 + cos(RotAz + deg2rad(240))*Blade(:,1);
    
    ublade_1 = powerLaw(HHwind,HubHt,zblade_1,alpha);
    ublade_2 = powerLaw(HHwind,HubHt,zblade_2,alpha);
    ublade_3 = powerLaw(HHwind,HubHt,zblade_3,alpha);
    
    wind_vector = [ublade_1,ublade_2,ublade_3];
end