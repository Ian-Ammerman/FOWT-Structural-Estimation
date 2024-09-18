% This is an example for plotting OpenFAST results

%% ----- Load OpenFAST Results ----- %%
load('OpenFAST/Simulations/5MW_Baseline_Fixed/Turbulent_13mps_Monhegan/5MW_Baseline_Fixed_FAST_Results.mat');

%% ----- Do Plotting ----- %%
figure
gca; hold on; box on;
plot(sim_results.Time, sim_results.TTDspFA)
title('Tower-Top Fore-Aft Displacement')
xlabel('Time [s]')
ylabel('Displacement [m]')
xlim([0,120])