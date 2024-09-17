clear all; close all; clc;

global btable ttable sim_results

%% ----- Load in Tower & Blade Information ----- %%
[btable, ttable] = load5MWDistributed;
blade = readmatrix('C:\Umaine Google Sync\Masters Working Folder\Multibody FOWT Model\Aero Data\NREL_5MW_Blade.csv');

%% ----- Load in OpenFAST Results ----- %%
% Steady 13 m/s simulation
% load("C:\Umaine Google Sync\Masters Working Folder\1 - OpenFAST\Simulations\5MW_Baseline_Fixed\SteadyWind_13mps\5MW_Baseline_Fixed_FAST_Results.mat");

% Turbulent 13 m/s simulation
load('C:\Umaine Google Sync\Masters Working Folder\1 - OpenFAST\Simulations\5MW_Baseline_Fixed\Turbulent_13mps_Monhegan\5MW_Baseline_Fixed_FAST_Results.mat')

% Unwrap rotor azimuth angle
rotAz = unwrap(sim_results.Azimuth);

%% ----- Run Standard Simulation ----- %%
% Initial conditions
x0 = zeros(8,1);
x = x0;

% Time vector
dt = 0.0125;
t = [0:dt:120];
N = length(t);

% Initialize output storage
x1store = zeros(length(x0), length(t));
x2store = zeros(length(x0), length(t));

% Applied forces structure
F = struct();

% Observer values
Q = 10^6*eye(8);
R = 10^-6*eye(1);

x1 = x;
x2 = x;

tic
flaplog = zeros(3,length(t));
for i = 1:length(t)
    % tic
    % Update azimuth measurement
    F.azimuth = deg2rad(rotAz(i));

    % Update rotor speed measurement
    F.omega = sim_results.RotSpeed(i) * (pi/30);
    RPM = sim_results.RotSpeed(i);

    % Get modal velocity
    dqb = [x(6), x(7), x(8)];

    % Get hub height wind speed
    HH_wind = sim_results.Wind1VelX(i) - x(5);

    % Project wind onto blades
    azimuth = deg2rad(sim_results.Azimuth(i));

    zblade_1 = 90 + cos(azimuth)*blade(:,1);
    zblade_2 = 90 + cos(azimuth + deg2rad(120))*blade(:,1);
    zblade_3 = 90 + cos(azimuth + deg2rad(240))*blade(:,1);

    ublade_1 = powerLaw(HH_wind,90,zblade_1,0.18);
    ublade_2 = powerLaw(HH_wind,90,zblade_2,0.18);
    ublade_3 = powerLaw(HH_wind,90,zblade_3,0.18);
    
    wind_vector = [ublade_1,ublade_2,ublade_3];

    % Get blade pitch commands
    blade_pitch = (pi/180)*[sim_results.BldPitch1(i), sim_results.BldPitch1(i), sim_results.BldPitch1(i)];

    % Get flapwise loads
    aero_forces = RT_BEM_Ning_v2(wind_vector,dqb,RPM,blade_pitch);

    % Applied forces
    F.forces = zeros(4,1);
    F.forces(1) = sum(aero_forces.Thrust);
    F.forces(2:4) = aero_forces.FlapwiseForce;

    flaplog(:,i) = aero_forces.FlapwiseForce;

    % F.forces = zeros(4,1);
    % F.forces(1) = sum(flapNet(i,:));
    % F.forces(2:4) = flapNet(i,:)';

    % Measurements for Kalman Filter
    measurements = [sim_results.QD_TFA1(i)];

    x1 = RK4(dt, x1, @FixedTrinityRHS, F);
    x2 = SRCUKF(dt, x2, @FixedTrinityRHS, @FixedTrinityMeasurement, Q, R, measurements, F);
    % x2 = EKF(dt, x2, @FixedTrinityRHS, @FixedTrinityJacobian, @FixedTrinityMeasurement, Q, R, measurements, F);

    % Store state vector
    x1store(:,i) = x1;
    x2store(:,i) = x2;

    % Store outputs
    y1store(:,i) = FixedTrinityMeasurement(x1,F);
    % y2store(:,i) = FixedTrinityMeasurement(x2,F);
    % toc
end
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
