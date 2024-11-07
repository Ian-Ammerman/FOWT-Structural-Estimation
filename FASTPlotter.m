% This is an example for plotting OpenFAST results
clear all; close all; clc;

%% ----- User Inputs ----- %%
model_name = '5MW_OC4Semi_WSt_WavesWN';
sim_name = 'Elsevier_LC03';

%% ----- Load OpenFAST Results ----- %%
load_path = sprintf('OpenFAST/Simulations/%s/%s/%s_FAST_Results.mat', model_name, sim_name, model_name);
load(load_path);

%% ----- Do Plotting ----- %%
% Results are stored in a structure called "sim_results", with field names
% corresponding to OpenFAST output names as they appear in the Outlist. THe
% first field is sim_results.Time which is simulation time in seconds.

% Plot tower fore-aft displacement
figure
gca; hold on; box on;
plot(sim_results.Time, sim_results.TTDspFA)
title('Tower-Top Fore-Aft Displacement')
xlabel('Time [s]')
ylabel('Displacement [m]')
xlim([0,120])