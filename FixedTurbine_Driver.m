clear all; close all; clc;

%% ----- Define Global Variables ----- %%
% btable and ttable are tables of constant blade property values
% (distributed mass, stiffness, etc.). Using global variables to easily
% pass through to RHS function. This could also be done through the F
% structure as well.
global btable ttable sim_results

%% ----- Load in Tower & Blade Information ----- %%
% Blade & tower structural properties
[btable, ttable] = load5MWDistributed;

% Blade aerodynamic values
blade = readmatrix('OpenFAST\5MW_Baseline\AeroData\NREL_5MW_Blade.csv');

%% ----- Load in OpenFAST Results ----- %%
% Turbulent 13 m/s simulation
load('OpenFAST\Simulations\5MW_Baseline_Fixed\Turbulent_13mps_Monhegan\5MW_Baseline_Fixed_FAST_Results.mat');

% Unwrap rotor azimuth angle
rotAz = unwrap(sim_results.Azimuth);

% Get tower velocity measurement
TTVelFA = gradient(sim_results.TTDspFA, sim_results.Time);
TTVelFA = lowpass(TTVelFA,1,1/mean(diff(sim_results.Time)));

% Lowpass filter UKF measurements
sim_results.RootMyb1 = lowpass(sim_results.RootMyb1,0.7,1/mean(diff(sim_results.Time)));
sim_results.RootMyb2 = lowpass(sim_results.RootMyb2,0.7,1/mean(diff(sim_results.Time)));
sim_results.RootMyb3 = lowpass(sim_results.RootMyb3,0.7,1/mean(diff(sim_results.Time)));

%% ----- Simulation Setup ----- %%
% Initial conditions
x0 = zeros(8,1);
x1 = x0;
x2 = x0;

% Time vector
dt = 0.01;
t = [0:dt:30];
N = length(t);

% Initialize output storage
x1store = zeros(length(x0), length(t)); % Baseline linear model
x2store = zeros(length(x0), length(t)); % Unscented Kalman filter

% Applied forces structure
F = struct();

% Observer values
Q = 10^6*eye(8);
R = 10^-6*eye(4);

%% ----- Run Simulation ----- %%
tic
for i = 1:length(t)
    % Update system measurements
    F.azimuth = deg2rad(rotAz(i));
    RPM = sim_results.RotSpeed(i);
    F.omega = RPM * (pi/30);
    
    % Measurements for Kalman Filter
    kf_measurements = [TTVelFA(i);
                       sim_results.RootMyb1(i)*10^-3;
                       sim_results.RootMyb2(i)*10^-3;
                       sim_results.RootMyb3(i)*10^-3];
    
    % Update BEM Inputs
    HH_wind = sim_results.Wind1VelX(i) - x2(5);
    wind_vector = projectWind(HH_wind, F.azimuth, blade, 90, 0.18);
    
    % Compute control inputs
    blade_pitch = (pi/180)*[sim_results.BldPitch1(i), sim_results.BldPitch1(i), sim_results.BldPitch1(i)];

    % Call BEM routine for aero forces
    aero_forces = RT_BEM_Ning_v2(wind_vector,RPM,blade_pitch);

    % Store RHS forces for solver
    F.forces = zeros(4,1);
    F.forces(1) = sum(aero_forces.Thrust);
    F.forces(2:4) = aero_forces.FlapwiseForce;

    % Perform state update
    x1 = RK4(dt, x1, @FixedTrinityRHS, F);
    x2 = SRCUKF(dt, x2, @FixedTrinityRHS, @FixedTrinityMeasurement, Q, R, kf_measurements, F);

    % Store state vector
    x1store(:,i) = x1;
    x2store(:,i) = x2;

    % Print elapsed time
    if rem(t(i),10) == 0
        fprintf('Elapsed Time: %0.3g seconds \n',t(i))
    end
end

% Compute simulation time & report efficiency
sim_time = toc;
fprintf('Elapsed Time: %0.3g | CPU Ratio: %0.3g \n',sim_time, sim_time/max(t));

%% Plot results
% Set plot limits
t2 = max(t);
t1 = max(t) - 30;

figure
gca; hold on; box on
plot(t, x1store(1,:), 'DisplayName', 'Linear', 'LineWidth', 1)
plot(t, x2store(1,:), 'DisplayName', 'SRCUKF', 'LineWidth', 1)
plot(sim_results.Time,sim_results.TTDspFA, 'LineStyle', '--', 'DisplayName','FAST', 'LineWidth', 1)
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Tower Fore-Aft Modal Displacement')
xlim([t1,t2])
legend

figure
gca; hold on; box on;
plot(t, x1store(2,:), 'DisplayName', 'Linear')
plot(t, x2store(2,:), 'DisplayName', 'SRCUKF')
plot(sim_results.Time,sim_results.Q_B1F1, 'DisplayName', 'FAST', 'LineStyle', '--', 'Color', 'Black')
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Blade 1 Flapwise Modal Displacement')
xlim([t1,t2])
legend

figure
gca; hold on; box on;
plot(t, gradient(x1store(5,:),t), 'DisplayName', 'Linear')
plot(t, gradient(x2store(5,:),t), 'DisplayName', 'SRCUKF')
plot(sim_results.Time,gradient(sim_results.QD_TFA1,sim_results.Time), 'DisplayName', 'FAST')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
title('Tower Fore-Aft Acceleration')
xlim([t1,t2])
legend

% figure
% gca; hold on; box on;
% yyaxis left
% plot(t, y1store(2,:), 'DisplayName', 'Linear')
% yyaxis right
% plot(sim_results.Time, sim_results.RootMyb1*10^3, 'DisplayName', 'FAST')
% xlabel('Time [s]')
% ylabel('Bending Moment [N-m]')
% title('Blade 1 Flapwise Root Bending Moment')
% xlim([min(t),max(t)])
% legend


% plot(t, xstore(3,:), 'DisplayName', 'Blade 2')
% plot(t, xstore(4,:), 'DisplayName', 'Blade 3')











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