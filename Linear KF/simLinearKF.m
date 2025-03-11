clear all; close all; clc;


% load('C:\Umaine Google Sync\GitHub\SHARK\OpenFAST\Simulations\FOCAL_C4\Steady_Wind_13ms\FOCAL_C4_FAST_Results.mat')
load('C:\Umaine Google Sync\GitHub\SHARK\OpenFAST\Simulations\FOCAL_C4\Elsevier_LC02\FOCAL_C4_FAST_Results.mat');

%% ----- Setup ----- %%
% Define constants
dt = 0.01;      % sample time [s]
Fs = 1/dt;      % sampling frequency
g = 9.806;      % gravitational acceleration
rt = 160.23;    % height from SWL of nacelle accelerometer

t1 = 14001;
t2 = 26001;


% Load FOCAL Data
t = trimData(0.01);

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
q = 10^0;
r = 10^-1;

Q = q*eye(4);
R = r*eye(3);

%% ----- Define System Measurements ----- %%
% Compute pitch acceleration from pitch velocity
t.PitchAcc = gradient(t.PitchVel*(pi/180),dt);

% Compute tower-top velocity
ttacc = rMean(t.accelNacelleAx + g*t.Pitch*(pi/180));
% ttvel = rMean(cumtrapz(t.Time,highpass(ttacc,0.05,1/dt,'ImpulseResponse','iir'))) - 0*rt*t.PitchVel*(pi/180);
ttvel = rMean(cumtrapz(t.Time,highpass(ttacc,0.15,1/dt,'ImpulseResponse','iir')));

% savetime = t.Time;
% save('Data/LC02_TTVel_Time.mat','savetime');
% save('Data/LC02_TTVel.mat','ttvel');

figure
gca; hold on; box on;
plot(sim_results.Time-700,gradient(sim_results.QD_TFA1,dt))
plot(t.Time,ttacc,'DisplayName','Accelerometer')
legend

figure
gca; hold on;
plot(t.Time,ttvel,'DisplayName','Integrated')
plot(sim_results.Time-700,sim_results.QD_TFA1,'DisplayName','OpenFAST')
legend
% plot(sim_results.Time-700,sim_results.NcIMUTVxs)

% Compute tower-top position from acceleration
ttpos = cumtrapz(t.Time,rMean(ttvel));

%% ----- Simulate Kalman Filter ----- %%
% Initial conditions
x0 = zeros(4,1);
xc = x0;

% Inputs are zero
u = zeros(2,1);

% Initialize output handling
kf_results = zeros(length(t.Time),length(x0));

% Initialize integrations
% ttvel = 0;
% ttpos = 0;
meas_vals = zeros(length(t.Time),3);
% Simulate
for i = 1:size(kf_results,1)

    % Get "gyro" measurements
    ddtheta_gyro = t.PitchAcc(i);
    dtheta_gyro = t.PitchVel(i) * (pi / 180);
    theta_gyro = t.Pitch(i) * (pi / 180);

    % % Get accelerometer measurement
    % a_raw = t.accelNacelleAx(i);
    % a_adj = a_raw + g*xc(1);
    % 
    % % Filter acceleration measurement
    % % [a_hp, highpass_state] = RTFilter(a_adj, highpass_coeffs, highpass_state);
    % 
    % % Correct for pitch acceleration
    % a_flx = a_adj + rt*ddtheta_gyro;
    % 
    % % Integrate acceleration
    % ttvel = ttvel + a_flx*dt;
    % 
    % % Integrate tower position
    % ttpos = ttpos + ttvel*dt;

    % Combine measurements
    z = [theta_gyro;dtheta_gyro;ttvel(i)];
    meas_vals(i,:) = z';

    % Inputs
    u(1) = rt*t.towerTopFx(i);
    u(2) = t.towerTopFx(i);

    % Compute corrected state vector
    xc = KF(dt,xc,u,z,Q,R,dsys);

    % Store results
    kf_results(i,:) = xc';

end

% KF = []
KF = [t.Time,kf_results];
% save('Data/KF_LC01.mat',"KF")

% %%
% figure
% gca; hold on; box on;
% plot(t.Time,kf_results(:,1)*(180/pi),'DisplayName','KF Estimate')
% plot(t.Time,t.Pitch,'DisplayName','Experiment')
% legend
% title('Platform Pitch Angle')

%%
fig = figure;
ax = gca; hold on; box on;
ax = setAxesInfo(ax);
xlim([160,200])

plot(sim_results.Time-700,sim_results.TTDspFA,'DisplayName','OpenFAST','Color','Black','LineWidth',1.5)
% plot(t.Time,ttpos,'DisplayName','Integrated')
plot(t.Time,kf_results(:,2),'DisplayName','KF Estimate','LineWidth',1.5,'LineStyle',':')
legend
xlabel('Time [s]')
ylabel({'Tower 1^{st} Mode Displacement [m]'})
% title('Floating Tower 1^{st} Mode Displacement')
% exportgraphics(fig,'Figures/Experimental_TTDspFA_TDComparison.pdf','ContentType','vector')
% 
% %%
% fig = figure;
% ax = gca; hold on; box on;
% ax = setAxesInfo(ax);
% xlim([520,550])
% 
% plot(t.Time,t.PitchVel*(pi/180),'DisplayName','Experiment','LineWidth',1.5,'Color','Black')
% plot(t.Time,kf_results(:,3),'DisplayName','KF Estimate','LineWidth',1.5,'LineStyle',':')
% 
% legend
% xlabel('Time [s]')
% ylabel('Platform Pitch Velocity [rads/s]')
% % title('Platform Pitch Velocity')
% % exportgraphics(fig,'Figures/Experimental_PitchVelocity_TDComparison.pdf','ContentType','vector')
% % % 
% %%
% fig = figure;
% ax = gca; hold on; box on;
% ax = setAxesInfo(ax);
% xlim([160,200])
% plot(t.Time,meas_vals(:,3),'DisplayName','Experiment','LineWidth',1.5,'Color','Black')
% plot(t.Time,kf_results(:,4),'DisplayName','KF Estimate','LineWidth',1.5,'LineStyle',':')
% % plot(t.Time,meas_vals(:,3),'DisplayName','Experiment','LineWidth',1.5,'Color','Black')
% legend
% ylabel('Tower-Top Velocity [m/s]')
% xlabel('Time [s]')
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
%% ----- Check PSD ----- %%
% epsd = myPSDv1(kf_results(:,2),Fs,0);
% fpsd = myPSDv1(sim_results.Q_TFA1,1/mean(diff(sim_results.Time)),0);

epsd = myFFT(rMean(kf_results(:,2)),Fs,0);
fpsd = myFFT(rMean(sim_results.Q_TFA1),1/mean(diff(sim_results.Time)),0);

fig = figure;
ax = gca; hold on; box on;
ax = setAxesInfo(ax);
set(gca, 'YScale', 'log')
xlim([0.01,0.25])
plot(fpsd(:,1),fpsd(:,2),'DisplayName','OpenFAST','Color','Black')
plot(epsd(:,1),epsd(:,2),'DisplayName','KF','LineWidth',1)

xlabel('Frequency [Hz]')
ylabel({'Floating Tower 1^{st} Mode','Displacement FFT [m]'})
% title({'Floating Tower 1^{st} Mode','Displacement FFT [m]'})
legend
% exportgraphics(fig,'Figures/Experimental_TTDspFA_FDComparison.pdf','ContentType','vector')