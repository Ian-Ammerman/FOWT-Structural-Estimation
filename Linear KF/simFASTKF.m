clear all; close all; clc;


load('C:\Umaine Google Sync\GitHub\SHARK\OpenFAST\Simulations\FOCAL_C4\Steady_Wind_13ms\FOCAL_C4_FAST_Results.mat')
load('C:\Umaine Google Sync\GitHub\SHARK\OpenFAST\Simulations\FOCAL_C4\Elsevier_LC02\FOCAL_C4_FAST_Results.mat');

%% ----- Setup ----- %%
% Define constants
dt = 0.05;      % sample time [s]
Fs = 1/dt;      % sampling frequency
g = 9.806;      % gravitational acceleration
rt = 145.23;    % height from SWL of nacelle accelerometer

t1 = 14001;
t2 = 26001;


% Load FOCAL Data
% % t = trimData(0.01);

%% ----- Define Dynamics Model ----- %%
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

% State-space model
A = [zeros(size(M)),eye(size(M));
     -M\K, -M\C];

B = [zeros(2);M\eye(2)];

% C = eye(4);
C = [1,0,0,0;
     0,0,1,0;
     0,0,0,1];

D = zeros(size(C,1),size(B,2));

sys = ss(A,B,C,D);
dsys = c2d(sys,dt);

%% ----- Define Kalman Filter Weights ----- %%
q = 10^-1;
r = 10^-2;

Q = q*diag([1,0.01,1,1]);
R = r*diag([1,1,1]);

%% ----- Define System Measurements ----- %%
ttvel = lowpass(sim_results.NcIMUTVxs - rt*sim_results.QD_P - sim_results.QD_Sg,0.3,'ImpulseResponse','iir');
pitch_vel = lowpass(sim_results.NcIMURVys*(pi/180),0.15,Fs,'ImpulseResponse','iir');

figure
gca; hold on; box on;
plot(sim_results.Time,sim_results.QD_P,'DisplayName','OpenFAST')
plot(sim_results.Time,pitch_vel,'DisplayName','Filtered IMU')
legend

% Check IMU measurement results
% figure
% gca; hold on; box on;
% plot(sim_results.Time,sim_results.QD_TFA1,'DisplayName','OpenFAST')
% % plot(sim_results.Time,sim_results.NcIMUTVxs,'DisplayName','Raw IMU')
% plot(sim_results.Time,ttvel,'DisplayName','Filtered')
% legend
% 
% f = myFFT(sim_results.NcIMUTVxs,Fs,0);
% f2 = myFFT(ttvel,Fs,0);
% 
% figure
% gca; hold on; box on;
% plot(f(:,1),f(:,2),'DisplayName','Raw IMU')
% plot(f2(:,1),f2(:,2),'DisplayName','Filtered IMU')
% 
% 
% % ttvel = sim_results.QD_TFA1 + rt*sim_results.QD_P;
% % pitch_vel = sim_results.QD_P;
% 
% % figure
% % gca; hold on; box on;
% % plot(sim_results.Time,highpass(sim_results.QD_TFA1,0.32,Fs,'ImpulseResponse','iir'),'DisplayName','OpenFAST')
% % % plot(sim_results.Time,sim_results.NcIMUTVxs,'DisplayName','IMU Raw')
% % plot(sim_results.Time,ttvel,'DisplayName','IMU Filtered')
% % legend
% 
% % figure
% % gca; hold on; box on;
% % plot(sim_results.Time,sim_results.QD_P,'DisplayName','OpenFAST')
% % plot(sim_results.Time,sim_results.NcIMURVys*(pi/180),'DisplayName','IMU Raw')
% % plot(sim_results.Time,pitch_vel*(pi/180),'DisplayName','IMU Filtered')
% % legend

%% ----- Simulate Kalman Filter ----- %%
% Initial conditions
x0 = zeros(4,1);
xc = x0;

% Inputs are zero
u = zeros(2,1);

% Initialize output handling
kf_results = zeros(length(sim_results.Time),length(x0));

% Initialize integrations
% ttvel = 0;
% ttpos = 0;
meas_vals = zeros(length(sim_results.Time),3);
% Simulate
for i = 1:size(kf_results,1)

    % Get "gyro" measurements
    dtheta_gyro = pitch_vel(i);% * (pi / 180);
    theta_gyro = sim_results.Q_P(i);

    % Combine measurements
    z = [theta_gyro;dtheta_gyro;ttvel(i)];
    meas_vals(i,:) = z';

    % Inputs
    u(1) = rt*sim_results.RotThrust(i)*10^3;
    u(2) = sim_results.RotThrust(i)*10^3;

    % Compute corrected state vector
    xc = KF(dt,xc,u,z,Q,R,dsys);

    % Store results
    kf_results(i,:) = xc';

end

% KF = []
% KF = [t.Time,kf_results];
% save('Data/KF_LC02.mat',"KF")

%%
figure
gca; hold on; box on;
plot(sim_results.Time,kf_results(:,1)*(180/pi),'DisplayName','KF Estimate')
plot(sim_results.Time,sim_results.PtfmPitch,'DisplayName','Experiment')
legend
title('Platform Pitch Angle')
% % 
%%
fig = figure;
ax = gca; hold on; box on;
ax = setAxesInfo(ax);
xlim([860,920])

plot(sim_results.Time,sim_results.Q_TFA1,'DisplayName','OpenFAST','Color','Black','LineWidth',1.5)
plot(sim_results.Time,kf_results(:,2),'DisplayName','KF Estimate','LineWidth',1.5,'LineStyle',':')
% plot(sim_results.Time,lowpass(kf_results(:,2),0.12,1/mean(diff(sim_results.Time)),'ImpulseResponse','iir'),'DisplayName','KF Filtered','LineWidth',1.5,'LineStyle','--')
legend
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Floating Tower 1^{st} Mode Displacement')
% exportgraphics(fig,'Figures/Experimental_TTDspFA_TDComparison.pdf','ContentType','vector')
% 
%%
fig = figure;
ax = gca; hold on; box on;
ax = setAxesInfo(ax);
xlim([520,550])

plot(sim_results.Time,sim_results.QD_P,'DisplayName','Experiment','LineWidth',1.5,'Color','Black')
plot(sim_results.Time,kf_results(:,3),'DisplayName','KF Estimate','LineWidth',1.5,'LineStyle',':')

legend
xlabel('Time [s]')
ylabel('Angular Velocity [rads/s]')
title('Platform Pitch Velocity')
% % % exportgraphics(fig,'Figures/Experimental_PitchVelocity_TDComparison.pdf','ContentType','vector')
% % 
%%
figure
gca; hold on; box on;
plot(sim_results.Time,kf_results(:,4),'DisplayName','KF Estimate')
plot(sim_results.Time,sim_results.QD_TFA1,'DisplayName','OpenFAST')
% plot(sim_results.Time,sim_results.NcIMUTVxs,'DisplayName','IMU Raw')
% plot(sim_results.Time,ttvel,'DisplayName','IMU Filtered')
legend
title('Tower-Top Velocity')
% 
% % figure('Position',[344.2,921,1094.4,420])
% % subplot(1,2,1)
% % gca; hold on; box on;
% % plot(t.Time,kf_results(:,2)*(180/pi),'DisplayName','KF Estimate')
% % plot(t.Time,t.Roll,'DisplayName','Experiment')
% % legend
% % title('Platform Roll Angle')
% % 
% % % figure
% % subplot(1,2,2)
% % gca; hold on; box on;
% % plot(t.Time,kf_results(:,4)*(180/pi),'DisplayName','KF Estimate')
% % plot(t.Time,t.RollVel,'DisplayName','Experiment')
% % legend
% % title('Platform Roll Velocity')
% 
% % %% ----- Check FFT ----- %%
% % fkf = myFFT(kf_results(:,2),1/dt);
% % fex = myFFT(sim_results.TTDspFA,1/mean(diff(sim_results.Time)));
% % 
% % figure
% % gca; hold on; box on;
% % plot(fkf(:,1),fkf(:,2),'DisplayName','KF')
% % plot(fex(:,1),fex(:,2),'DisplayName','FAST')
% % xlim([0,1])
% 
% %% ----- Check PSD ----- %%
% % epsd = myPSDv1(kf_results(:,2),Fs,0);
% % fpsd = myPSDv1(sim_results.Q_TFA1,1/mean(diff(sim_results.Time)),0);
% 
% epsd = myFFT(rMean(kf_results(:,2)),Fs,0);
% fpsd = myFFT(rMean(sim_results.Q_TFA1),1/mean(diff(sim_results.Time)),0);
% 
% fig = figure;
% ax = gca; hold on; box on;
% ax = setAxesInfo(ax);
% set(gca, 'YScale', 'log')
% xlim([0.01,0.25])
% plot(fpsd(:,1),fpsd(:,2),'DisplayName','OpenFAST','Color','Black')
% plot(epsd(:,1),epsd(:,2),'DisplayName','KF','LineWidth',1.5)
% 
% xlabel('Frequency [Hz]')
% ylabel('Displacement Amplitude [m]')
% title('Floating Tower 1^{st} Mode Displacement FFT')
% legend
% exportgraphics(fig,'Figures/Experimental_TTDspFA_TDComparison.pdf','ContentType','vector')