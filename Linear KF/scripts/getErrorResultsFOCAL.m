clear all; close all; clc;

%% ===== Load Data ===== %%
% LC01
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\FOCAL Pitch Estimation\Data\FOCAL_4.2_M02D01T03_AC01_E02W01_R01_Z1_A1.mat')
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\FOCAL Pitch Estimation\Data\KF_LC01.mat')
% load('C:\Umaine Google Sync\GitHub\SHARK\OpenFAST\Simulations\FOCAL_C4\Elsevier_LC01\FOCAL_C4_FAST_Results.mat');
% % load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\FOCAL Pitch Estimation\Data\FOCAL_4.2_M02D01T01_AC01_E13W02_R01_Z1_A1.mat')
% t = trimData(mean(diff(sim_results.Time)));

% LC02
load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\FOCAL Pitch Estimation\Data\KF_LC02.mat')
load('C:\Umaine Google Sync\GitHub\SHARK\OpenFAST\Simulations\FOCAL_C4\Elsevier_LC02\FOCAL_C4_FAST_Results.mat');
load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\FOCAL Pitch Estimation\Data\FOCAL_4.2_M02D01T01_AC01_E13W02_R01_Z1_A1.mat')
t = trimData(mean(diff(sim_results.Time)));

% % LC03
% load('C:\Umaine Google Sync\GitHub\SHARK\OpenFAST\Simulations\FOCAL_C4\Elsevier_LC03\FOCAL_C4_FAST_Results.mat');
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\FOCAL Pitch Estimation\Data\FOCAL_4.2_M02D01T05_AC01_E23W03_R01_Z1_A1.mat')
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\FOCAL Pitch Estimation\Data\KF_LC03.mat')
% t = trimData(mean(diff(sim_results.Time)));

figure
gca; hold on; box on;
plot(t.Time,t.TTDspFA,'DisplayName','Raw');
plot(t.Time,lowpass(t.TTDspFA,0.3,1/mean(diff(t.Time)),'ImpulseResponse','iir'),'DisplayName','Filtered');

%% ===== Prepare Baseline ===== %%
dt = mean(diff(sim_results.Time));
g = 9.806;
rt = 160.23; 
% Compute pitch acceleration from pitch velocity
t.PitchAcc = gradient(t.PitchVel*(pi/180),dt);

% Compute tower-top velocity
ttacc = rMean(t.accelNacelleAx + g*t.Pitch*(pi/180));
% ttvel = rMean(cumtrapz(t.Time,highpass(ttacc,0.05,1/dt,'ImpulseResponse','iir'))) - 0*rt*t.PitchVel*(pi/180);
ttvel = rMean(cumtrapz(t.Time,highpass(ttacc,0.15,1/dt,'ImpulseResponse','iir')));
ttpos = cumtrapz(t.Time,ttvel);

% Trim time to 700 seconds
idx1 = dsearchn(channels(:,1),700);
idx2 = dsearchn(channels(:,1),1300);
channels = channels(idx1:idx2,:);
% ttvel = ttvel(idx1:idx2);
% savetime = savetime(idx1:idx2)-savetime(idx1);
channels(:,1) = channels(:,1) - min(channels(:,1));

idx1 = dsearchn(sim_results.Time,700);
idx2 = dsearchn(sim_results.Time,1300);
fields = fieldnames(sim_results);
for i = 1:length(fields)
    sim_results.(fields{i}) = sim_results.(fields{i})(idx1:idx2);
end
sim_results.Time = sim_results.Time - min(sim_results.Time);

% Interpolate experiment to match FAST
newchannels = zeros(length(sim_results.Time),size(channels,2));
newchannels(:,1) = sim_results.Time;
for i = 2:size(channels,2)
    newchannels(:,i) = pchip(channels(:,1),channels(:,i),sim_results.Time);
end

% Interpolate kalman filter
new_kf = zeros(length(sim_results.Time),4);
for i = 1:4
    new_kf(:,i) = pchip(KF(:,1),KF(:,i+1),sim_results.Time);
end

% figure
% gca; hold on;
% plot(sim_results.Time,sim_results.PtfmPitch,'DisplayName','FAST');
% plot(newchannels(:,1),newchannels(:,9),'DisplayName','Interp. Data');
% plot(channels(:,1),channels(:,9),'DisplayName','Raw Data')
% legend

% Full baseline
baseline = [t.Pitch*(pi/180),lowpass(t.TTDspFA,0.3,1/mean(diff(t.Time)),'ImpulseResponse','iir'),t.PitchVel*(pi/180),ttvel];

figure
gca; hold on; box on;
plot(sim_results.Time,sim_results.TTDspFA,'DisplayName','OpenFAST')
plot(sim_results.Time,baseline(:,2),'DisplayName','From 6DOF')
plot(sim_results.Time,new_kf(:,2),'DisplayName','KF')
title('Tower Displacement')
legend
% 
% figure
% gca; hold on; box on;
% plot(sim_results.Time,baseline(:,4),'DisplayName','Baseline')
% plot(sim_results.Time,sim_results.QD_TFA1,'DisplayName','OpenFAST')
% plot(sim_results.Time,new_kf(:,4),'DisplayName','KF')
% % plot(sim_results.Time,gradient(new_kf(:,2),mean(diff(sim_results.Time))),'DisplayName','KF Gradient')
% title('Tower Velocity')
% legend
% 
% figure
% gca; hold on; box on;
% plot(sim_results.Time,highpass(sim_results.Q_TFA1,0.15,1/mean(diff(sim_results.Time))),'DisplayName','OpenFAST')
% plot(t.Time,baseline(:,2),'DisplayName','From 6DOF')
% title('Tower Position')
% legend

%% ===== Mean Relative Error ===== %%
all_error = 100*mre(new_kf,baseline)';
% save('Data/LC03_Error.mat',"all_error");

%% ===== Cross-Correlation ===== %%

all_R2 = diag(corr(new_kf,baseline)).^2;


% signal1 = baseline(:,2);
% signal2 = new_kf(:,2);
% % [xcorr_vals, lags] = xcorr(signal1, signal2, 'coeff');
% % [~, max_idx] = max(xcorr_vals);
% % phase_shift = lags(max_idx);
% % stem(xcorr_vals,lags);
% 
% r2 = getCorrelationCoefficient(signal1,signal2)