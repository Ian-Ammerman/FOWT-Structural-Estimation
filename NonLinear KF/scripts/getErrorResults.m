clear all; close all; clc;

% % LC01
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\Elsevier_LC01\5MW_OC4Semi_WSt_WavesWN_FAST_Results.mat');
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\EKF_LC01.mat');
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\EKF_wFBG_LC01.mat');
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\UKF_LC01.mat');
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\UKF_wFBG_LC01.mat');

% LC02
load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\Elsevier_LC02\5MW_OC4Semi_WSt_WavesWN_FAST_Results.mat');
load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\EKF_LC02.mat');
load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\EKF_wFBG_LC02.mat');
load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\UKF_LC02.mat');
load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\UKF_wFBG_LC02.mat');

% LC03
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\Elsevier_LC03\5MW_OC4Semi_WSt_WavesWN_FAST_Results.mat');
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\EKF_LC03.mat');
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\EKF_wFBG_LC03.mat');
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\UKF_LC03.mat');
% load('C:\Umaine Google Sync\Masters Working Folder\State Observation\Examples\NREL5MW_OC4_State_Estimator\Data\UKF_wFBG_LC03.mat');

% Isolate channels
baseline = [sim_results.Q_P,sim_results.Q_TFA1,sim_results.Q_B1F1,sim_results.Q_B2F1,sim_results.Q_B3F1,sim_results.QD_P,sim_results.QD_TFA1,sim_results.QD_B1F1,sim_results.QD_B2F1,sim_results.QD_B3F1];

%% ===== Mean Relative Error ===== %%
LC01_EKF_mre = 100*mre(EKF(:,2:end),baseline)';
LC01_EKFwFBG_mre = 100*mre(EKF_wFBG(:,2:end),baseline)';
LC01_UKF_mre = 100*mre(UKF(:,2:end),baseline)';
LC01_UKFwFBG_mre = 100*mre(UKF_wFBG(:,2:end),baseline)';

all_mre = [LC01_EKF_mre,LC01_EKFwFBG_mre,LC01_UKF_mre,LC01_UKFwFBG_mre];

% writematrix(all_mre,'Data/LC03_MRE.dat','FileType','text')

%% ===== Coefficient of Determination ===== %%
EKF_r2 = diag(corr(EKF(:,2:end),baseline)).^2;
EKFwFBG_r2 = diag(corr(EKF_wFBG(:,2:end),baseline)).^2;
UKF_r2 = diag(corr(UKF(:,2:end),baseline)).^2;
UKFwFBG_r2 = diag(corr(UKF_wFBG(:,2:end),baseline)).^2;

all_r2 = [EKF_r2,EKFwFBG_r2,UKF_r2,UKFwFBG_r2];

figure
gca; hold on; box on;
plot(sim_results.Time,sim_results.QD_B1F1,'DisplayName','FAST')
plot(UKF_wFBG(:,1),UKF_wFBG(:,9),'DisplayName','UKF');
plot(EKF_wFBG(:,1),EKF_wFBG(:,9),'DisplayName','EKF');
xlim([160,180])
legend