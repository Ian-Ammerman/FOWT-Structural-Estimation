clear all; close all; clc;

%% ----- Load in & Process Data ----- %%
% Load data
load('C:\Umaine Google Sync\GitHub\SHARK\OpenFAST\Simulations\FOCAL_C4\Elsevier_LC02\FOCAL_C4_FAST_Results.mat');

% Sampling info
t = sim_results.Time;
dt = mean(diff(t));
Fs = 1/dt;


%% ===== Tower Linear Velocity ===== %%
% Tower natural frequency
ftow = 0.385;
bprange = [0.32,0.5]

% Extract IMU measurement
tower_vel = highpass(sim_results.QD_TFA1,0.32,Fs,'ImpulseResponse','iir');
tower_imu_raw = sim_results.NcIMUTVxs - 87.6*sim_results.QD_P;
tower_imu_filtered = highpass(tower_imu_raw,0.32,Fs,'ImpulseResponse','iir');

% Integrate for tower deflection and compare
tower_pos_integrated = cumtrapz(sim_results.Time,tower_imu_filtered);
tower_pos_baseline = highpass(sim_results.Q_TFA1,0.32,Fs,'ImpulseResponse','iir');

% FFT of tower motions
ft = myFFT(tower_vel,Fs,5);
ftimu = myFFT(tower_imu_raw,Fs,5);
ftimuf = myFFT(tower_imu_filtered,Fs,5);
fpos = myFFT(tower_pos_integrated,Fs,5);

figure
plot(fpos(:,1),fpos(:,2))


% % Plot results
% figure
% subplot(1,2,1)
% gca; hold on; box on;
% plot(ft(:,1),ft(:,2),'DisplayName','Baseline')
% % plot(ftimu(:,1),ftimu(:,2),'DisplayName','IMU Raw')
% plot(ftimuf(:,1),ftimuf(:,2),'DisplayName','IMU Filtered')
% title('Tower Velocity FFT')
% xlim([0,0.6])
% legend
% 
% subplot(1,2,2)
% gca; hold on; box on;
% plot(sim_results.Time,tower_vel,'DisplayName','Baseline');
% % plot(sim_results.Time,tower_imu_raw,'DisplayName','IMU Raw');
% plot(sim_results.Time,tower_imu_filtered,'DisplayName','IMU Filtered')
% title('Tower Velocity Time-Series')
% legend

% Plot results
figure
gca; hold on; box on;
plot(sim_results.Time,tower_pos_baseline,'DisplayName','Baseline');
% plot(sim_results.Time,tower_imu_raw,'DisplayName','IMU Raw');
plot(sim_results.Time,rMean(tower_pos_integrated),'DisplayName','Integrated')
title('Tower Position Time-Series')
legend


% %% ===== Platform Pitch Velocity ===== %%
% % Extract IMU angular velocity measurement
% pitch_vel = sim_results.QD_P;
% imu_raw = sim_results.NcIMURVys*(pi/180);
% 
% fp = myFFT(pitch_vel,Fs,0);
% fimu = myFFT(imu_raw,Fs,0);
% 
% % Plot results
% figure
% subplot(1,2,1)
% gca; hold on; box on;
% plot(fp(:,1),fp(:,2),'DisplayName','Baseline')
% plot(fimu(:,1),fimu(:,2),'DisplayName','IMU')
% title('Platform Pitch Velocity')
% legend
% xlim([0,0.8])
% 
% subplot(1,2,2)
% gca; hold on; box on;
% plot(sim_results.Time,pitch_vel,'DisplayName','Baseline')
% plot(sim_results.Time,imu_raw,'DisplayName','Raw IMU')
% title('Platform Pitch Velocity')
% legend