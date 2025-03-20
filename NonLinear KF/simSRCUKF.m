clear all; close all; clc;

%% ===== Load in Simulation Data ===== %%
% Define simulation
sim_name = 'Elsevier_LC02';
load(sprintf('Data/%s/5MW_OC4Semi_WSt_WavesWN_FAST_Results.mat',sim_name));

% Simulation parameters
dt = mean(diff(sim_results.Time));      % sample time [s]
Fs = 1/dt;                              % sampling frequency
g = 9.806;                              % gravitational acceleration
rt = 90;                            % height from SWL of nacelle accelerometer

%% ===== Design Kalman Filter ===== %%
% Tuning matrix weights
q = 10^-2;
r = 10^-0;

Q = q*diag([1,0.1,100,100,100,  1,0.1,0.001,0.001,0.001]);
R = r*diag([1,1000,0.001]);
% R = r*diag([1,100,0.001,0.1,0.1,0.1]);

% Floating system properties
rotor = loadNREL5MWRotor;
tower = loadNREL5MWTower;
floater = loadNREL5MWFloater; 

% Pre-process measurement signals
sim_results.NcIMUTVxs = highpass(sim_results.NcIMUTVxs,0.2,1/dt,"ImpulseResponse",'iir');
sim_results.NcIMURVys = lowpass(sim_results.NcIMURVys,0.2,1/dt,'ImpulseResponse','iir');

%% ===== Run Simulation ===== %%
% Initial conditions
x0 = zeros(10,1);
xc = x0;

% Inputs are zero
u = zeros(2,1);

% Initialize output handling
kf_results = zeros(length(sim_results.Time),length(x0));
% meas_vals = zeros(length(sim_results.Time),3);
meas_vals = zeros(length(sim_results.Time),length(diag(R)));

% Simulate
for i = 1:size(kf_results,1)

    % Get rotor info
    azimuth = sim_results.Azimuth(i)*(pi/180) + [0,2*pi/3,4*pi/3];
    rotor.azimuth = azimuth;
    omega = sim_results.RotSpeed(i)*(pi/30);

    % Get "gyro" measurements
    ttvel = sim_results.NcIMUTVxs(i);
    dtheta_gyro = sim_results.NcIMURVys(i)*(pi/180);
    theta_gyro = sim_results.Q_P(i);
    
    % Combine measurements
    % z = [theta_gyro;dtheta_gyro;ttvel;sim_results.Q_B1F1(i);sim_results.Q_B2F1(i);sim_results.Q_B3F1(i)];
    % z = [theta_gyro;dtheta_gyro;ttvel;sim_results.QD_B1F1(i);sim_results.QD_B2F1(i);sim_results.QD_B3F1(i)];
    z = [theta_gyro;dtheta_gyro;ttvel];
    meas_vals(i,:) = z';

    % Compute blade loads
    Uinf = sim_results.Wind1VelX(i);    % free-stream wind velocity
    Urel = Uinf - xc(7) - tower.HH*xc(6);     % hub-height apparent wind speed
    wind_vector = projectWind(Urel, azimuth(1), rotor.blade, tower.HH, 0.18);
    rotor.RPM = sim_results.RotSpeed(i);
    rotor.dqb = [xc(8),xc(9),xc(10)];
    rotor.blade_pitch = (pi/180)*[sim_results.BldPitch1(i), sim_results.BldPitch2(i), sim_results.BldPitch3(i)];
    aero_forces = RT_BEM_Ning_v3(wind_vector,rotor);

    % Inputs
    u = zeros(5,1);
    u(1) = rt*sim_results.RotThrust(i)*10^3;
    u(2) = sim_results.RotThrust(i)*10^3;
    u(3:5) = aero_forces.FlapwiseForce;
    % u(3) = (1/3)*sim_results.RotThrust(i)*10^3*0.35;
    % u(4) = (1/3)*sim_results.RotThrust(i)*10^3*0.35;
    % u(5) = (1/3)*sim_results.RotThrust(i)*10^3*0.35;

    % Compute corrected state vector
    M = NREL5MW_MassMatrix(azimuth,tower,rotor,floater);
    K = NREL5MW_StiffnessMatrix(azimuth,omega,tower,rotor,floater);
    C = NREL5MW_DampingMatrix(K,M);

    RHSFunc = @(x,u) NREL5MW_RHSFunc(x,u,M,C,K);

    xc = SRCUKF(dt,xc,RHSFunc,@NREL5MW_MeasFunc,Q,R,z,u,0.5);

    % Store results
    kf_results(i,:) = xc';

end

%% ===== Process Results ===== %%
% Save results
UKF_wFBG = [sim_results.Time,kf_results];
save('Data/UKF_wFBG_LC02.mat','UKF_wFBG');

% UKF = [sim_results.Time,kf_results];
% save('Data/UKF_LC02.mat','UKF');

%% ===== Plot Results ===== %%
close all;

%%% Platform pitch angle
figure
ax = gca; hold on; box on;
ax = setAxesInfo(ax);

plot(sim_results.Time,kf_results(:,1)*(180/pi),'DisplayName','KF Estimate')
plot(sim_results.Time,sim_results.PtfmPitch,'DisplayName','Experiment')

legend
xlabel('Time [s]')
ylabel('Pitch Angle [deg]')
title('Platform Pitch Angle')

%%% Tower-Top Displacement
fig = figure;
ax = gca; hold on; box on;
ax = setAxesInfo(ax);
% xlim([190,240])

plot(sim_results.Time,sim_results.TTDspFA,'DisplayName','OpenFAST','LineWidth',1.5,'Color','Black')
plot(sim_results.Time,kf_results(:,2),'DisplayName','KF Estimate','LineWidth',1.5,'LineStyle',':')

legend
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Tower-Top Displacement')

%%% Blade 1 Flapwise Displacement
figure
gca; hold on; box on;

plot(sim_results.Time,sim_results.Q_B1F1,'DisplayName','OpenFAST','LineWidth',1.5,'Color','Black')
plot(sim_results.Time,kf_results(:,3),'DisplayName','KF Estimate','LineWidth',1.5,'LineStyle',':')

legend
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Blade 1 Flapwise Displacement')

%%% Platform Pitch Angular Velocity
figure
ax = gca; hold on; box on;
ax = setAxesInfo(ax);

plot(sim_results.Time,sim_results.QD_P,'DisplayName','OpenFAST','LineWidth',1.5,'Color','Black')
plot(sim_results.Time,kf_results(:,6),'DisplayName','EKF','LineWidth',1.5,'LineStyle',':')

legend
xlabel('Time [s]')
ylabel('Angular Velocity [rads/s]')
title('Platform Pitch Velocity')

%%% Tower-Top Linear Velocity (Mode 1)
figure
gca; hold on; box on;
plot(sim_results.Time,kf_results(:,7),'DisplayName','KF Estimate')
plot(sim_results.Time,sim_results.QD_TFA1,'DisplayName','Experiment')
legend
xlabel('Time [s]')
ylabel('Velocity [m/s]')
title('Tower-Top Velocity')

%%% Blade 1 Flapwise Velocity (Mode 1)
figure
gca; hold on; box on;
plot(sim_results.Time,kf_results(:,8),'DisplayName','KF Estimate')
plot(sim_results.Time,sim_results.QD_B1F1,'DisplayName','OpenFAST')
legend
xlabel('Time [s]')
ylabel('Velocity [m/s]')
title('Blade 1 Flapwise Velocity')


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