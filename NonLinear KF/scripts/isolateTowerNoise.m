clear all; close all; clc;

%% ----- Load in & Process Data ----- %%
% Load OpenFAST data
load('Data/Elsevier_LC02/5MW_OC4Semi_WSt_WavesWN_FAST_Results.mat');

% Sampling info
t = sim_results.Time;
dt = mean(diff(t));
Fs = 1/dt;

% Extract IMU angular velocity measurement
pitch_vel = sim_results.QD_P;
imu_raw = sim_results.NcIMURVys*(pi/180);

fp = myFFT(pitch_vel,Fs);
fimu = myFFT(imu_raw,Fs);

figure
gca; hold on; box on;
plot(fp(:,1),fp(:,2),'DisplayName','Baseline')
plot(fimu(:,1),fimu(:,2),'DisplayName','IMU')
xlim([0,0.8])





% 
% figure
% gca; hold on; box on;
% plot(sim_results.Time,pitch_vel,'DisplayName','Baseline')
% plot(sim_results.Time,imu_raw,'DisplayName','Raw IMU')