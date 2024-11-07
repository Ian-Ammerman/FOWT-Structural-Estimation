clear all; close all; clc;

%% ----- User Inputs ----- %%
% Simulation setup
x0 = zeros(12,1);               % initial conditions
x = x0;                         % initialize state vector
xdot = zeros(length(x0),1);     % Initialize state derivative

dt = 0.01;          % simulation time step (must match FAST)
tmax = 50;         % simulation stop time (must be <= TMax from FAST

% Observer setup
Q = 10^6*eye(length(x));    % Process uncertainty (for Kalman filter)
R = 10^-6*eye(5);           % Measurement uncertainty (for Kalman filter)

% "True" system for comparison
turbine_name = '5MW_Baseline';              % Turbine model, used to locate blade/tower props
system_name = '5MW_OC4Semi_WSt_WavesWN';    % Name of OpenFAST model used for comparison
sim_name = 'Elsevier_LC02';                 % Name of simulation to run

load(sprintf('OpenFAST/Simulations/%s/%s/%s_FAST_Results.mat', system_name, sim_name, system_name));

% System measurements
sys_measurements = {'Q_P', 'TwrBsMyt','RootMyb1','RootMyb2','RootMyb3'};     % FAST output channels to treat as "measurements"

for i = 1:length(sys_measurements)
    Y(:,i) = sim_results.(sys_measurements{i});
end

% Scale measuremnts
% Y(:,1) = Y(:,1)*(pi/180);
Y(:,2) = Y(:,2)*10^-5;
Y(:,3:5) = Y(:,3:5)*10^-3;

% % Simulated measurement noise (optional)
% add_noise = false;      % flag whether to add white noise to measurement channels
% noiseStd = [1,1];       % standard deviation of noise to add to each channel, in order
% seed = [1,1];           % seed for RNG within addWhiteNoise function
% 
% if add_noise
%     for i = 1:size(Y,2)
%         sig = Y(:,i);
%         Y(:,i) = addWhiteNoise(sig, noiseStd(i), seed(i));
%     end
% end

% Load "Ground-Truth Vectors
states = {'Q_P','Q_TFA1','Q_B1F1','Q_B2F1','Q_B3F1','Time','QD_P','QD_TFA1','QD_B1F1','QD_B2F1','QD_B3F1','Time'};
truth = zeros(length(sim_results.Time),length(states));
for i = 1:length(states)
    truth(:,i) = sim_results.(states{i});
end

%% ----- Initialize Floating System ----- %%
% Rotor information object
Rotor = loadNREL5MWRotor;   % load rotor information object for NREL 5MW reference turbine
% Rotor = loadFOCALRotor;     % load rotor information object for FOCAL model of IEA-15MW
% Rotor = loadIEA15MWRotor;   % load rotor information object for IEA-15MW reference turbine

% Tower information object
Tower = loadNREL5MWTower;   % load tower information object for NREL 5MW reference turbine
% Tower = loadFOCALTower;     % load tower information object for FOCAL model of IEA-15MW
% Tower = loadIEA15MWTower;   % load tower information object for IEA-15MW reference turbine

%% ----- Initialize Simulation ----- %%
% Simulation time vector
t = [0:dt:tmax];     % simulation time vector
N = length(t);      % length of simulation time vector (useful to have)

% Storage Variables
xstore = zeros(length(x0), length(t));

% Input structure
F = struct();
F.tower = Tower;    % pass tower object to inner functions
F.rotor = Rotor;    % pass rotor object to inner functions

% Define RHS function 
RHSFunc = @(x, F) RHS(x, F, @NREL5MW_MassMatrix, @NREL5MW_DampingMatrix, @NREL5MW_StiffnessMatrix);

NEES = zeros(length(t),1);
%% ----- Run Standard Simulation ----- %%
tic
for i = 1:length(t)

    % Environmental
    Uinf = sim_results.Wind1VelX(i);    % free-stream wind velocity
    Urel = Uinf - x(8) - 87.6*x(7);     % hub-height apparent wind speed

    % System "measurements"
    Rotor.RPM = sim_results.RotSpeed(i);            % Rotor speed [RPM]
    F.omega = Rotor.RPM * (pi/30);                  % Rotor speed [rads/s]
    F.azimuth = deg2rad(sim_results.Azimuth(i));    % Rotor azimuth [rads]
    Rotor.dqb = [x(9),x(10),x(11)];                       % Blade local velocities [m/s]
    
    % Measurements for Kalman Filter
    kf_measurements = Y(i,:)';
    
    % Project wind onto blades
    wind_vector = projectWind(Urel, F.azimuth, Rotor.blade, 90, 0.18);

    % Get blade pitch commands
    Rotor.blade_pitch = (pi/180)*[sim_results.BldPitch1(i), sim_results.BldPitch2(i), sim_results.BldPitch3(i)];

    % Get loads
    % aero_forces = RT_BEM_Ning_v2(wind_vector,dqb,RPM,blade_pitch);
    aero_forces = RT_BEM_Ning_v3(wind_vector,Rotor);
    % Ftwr_aero = TwrAeroRT(Uinf,x(8),0.18);
    % hydro_forces = RT_Hydro(dt, x, xdot);

    % Applied forces
    F.forces = zeros(6,1);
    F.forces(1) = sum(aero_forces.Thrust)*87.6;% + hydro_forces.PtfmPitch;
    F.forces(2) = sum(aero_forces.Thrust);% + Ftwr_aero;
    F.forces(3:5) = aero_forces.FlapwiseForce;
    F.forces(6) = 0;
    % F.forces(3:5)

    % Update state vector
    % xp = RK4(dt, x, RHSFunc, F);
    [xp,P] = SRCUKF(dt, x, RHSFunc, @NREL5MW_Measurement, Q, R, kf_measurements, F);

    % Compute NEES
    % NEES(i) = 

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
fprintf('Simulated Time: %0.3g \n',max(t)-min(t))
fprintf('CPU Time: %0.3g \n', tottime);
fprintf('CPU Time Ratio: %0.3g \n', tottime/(max(t)-min(t)));

%% ----- Plot Results ----- %%
close all
tmin = 0;
tmax = 100;

figure
ax = gca; hold on; box on;
plot(t, xstore(1,:),'DisplayName', 'Estimator','LineWidth',1)
plot(sim_results.Time,sim_results.PtfmPitch*(pi/180), 'DisplayName', 'FAST','Color','Black','LineStyle','--','LineWidth',1);
xlabel('Time [s]')
xlim([tmin,tmax])
ylabel('Displacement [rads]')
title('Platform Pitch Displacement')
legend

figure
gca; hold on; box on;
plot(t, xstore(2,:),'DisplayName', 'Estimator','LineWidth',1)
plot(sim_results.Time,sim_results.TTDspFA, 'DisplayName', 'FAST','Color','Black','LineStyle','--','LineWidth',1);
xlabel('Time [s]')
xlim([tmin,tmax])
ylabel('Displacement [m]')
title('Tower Top Displacement')
legend

figure
gca; hold on; box on;
plot(t, xstore(3,:),'DisplayName', 'Estimator','LineWidth',1)
plot(sim_results.Time,sim_results.Q_B1F1, 'DisplayName', 'FAST','Color','Black','LineStyle','--','LineWidth',1);
xlabel('Time [s]')
xlim([tmin,tmax])
ylabel('Displacement [m]')
title('Blade 1 Flapwise Modal Displacement')
legend

figure
gca; hold on; box on;
plot(t,xstore(7,:));
plot(sim_results.Time,truth(:,7))

figure
gca; hold on; box on;
plot(t,xstore(8,:));
plot(sim_results.Time,truth(:,8))

figure
gca; hold on; box on;
plot(t,xstore(9,:));
plot(sim_results.Time,truth(:,9))

% %% ----- Plot All References ----- %%
% for i = 1:length(states)
%     figure
%     gca; hold on; box on;
%     plot(t,xstore(i,:),'DisplayName','Estimate')
%     plot(t,truth(:,i),'DisplayName','OpenFAST')
%     xlabel('Time [s]')
%     title(states{i})
%     legend
% end

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