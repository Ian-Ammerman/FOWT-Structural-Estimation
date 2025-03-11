clear all; close all; clc;

%% ----- Load Inputs from OpenFAST ----- %%
% Load sim results
load('C:\Umaine Google Sync\Masters Working Folder\1 - OpenFAST\Simulations\FOCAL_C4\Steady_13mps_NoWave\FOCAL_C4_FAST_Results.mat');

% Extract useful quantities
time = sim_results.Time;
thrust = sim_results.RotThrust*10^3;
pitchm = sim_results.RotThrust*10^3*160.23;

% Simulation time step
dt = mean(diff(time));

%% ----- Perform Simulation ----- %%
q = zeros(4,1);
qstore = zeros(length(time),4);

for i = 1:length(sim_results.Time)
    q = RK4(dt,q,@rhs,[pitchm(i);thrust(i)]);
    qstore(i,:) = q';
end

figure
gca; hold on; box on;
plot(sim_results.Time,sim_results.Q_P,'DisplayName','OpenFAST')
plot(time,qstore(:,1),'DisplayName','Linear')
xlabel('Time [s]')
ylabel('Platform Pitch [rads]')
legend

figure
gca; hold on; box on;
plot(sim_results.Time,sim_results.Q_TFA1,'DisplayName','OpenFAST')
plot(time,qstore(:,2),'DisplayName','Linear')
xlabel('Time [s]')
ylabel('Tower Displacement [m]')
legend

figure
gca; hold on; box on;
plot(sim_results.Time,sim_results.QD_P,'DisplayName','OpenFAST')
plot(time,qstore(:,3),'DisplayName','Linear')
xlabel('Time [s]')
ylabel('Platform Pitch Rate [rads/s]')
legend

figure
gca; hold on; box on;
plot(sim_results.Time,sim_results.QD_TFA1,'DisplayName','OpenFAST')
plot(time,qstore(:,4),'DisplayName','Linear')
xlabel('Time [s]')
ylabel('Tower Velocity [m/s]')
legend

function xddot = rhs(x,u)
    % Mass matrix
    m11 = 48410870305.1658;
    m22 = 1278084.80047176;
    m12 = 185459490.208370;
    m21 = 185459490.208370;
    M = [m11, m12;
         m21, m22];
    
    % Damping matrix
    c1 = 20*197713958.925553;
    c2 = 10*20701.7415674758;
    C = [c1,  0;
          0, c2];
    
    % Stiffness matrix
    k1 = 4.0848e+09;  % Hydrostatic + mooring stiffness
    k2 = 3.8235e+06;  % Tower modal stiffness
    K = [k1,  0;
          0, k2];

    % Compute acceleration
    pos = x(1:size(x,1)/2,1);
    vel = x(size(x,1)/2+1:end,1);
    acc = M\(u - C*vel - K*pos);

    % Output
    xddot = [vel;acc];
    
end