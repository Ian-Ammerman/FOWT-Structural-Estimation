clear all; close all; clc;

%% ----- Load in Data ----- %%
% Load experimental results
load('../FOCAL Pitch Estimation/FOCAL_4.2_M02D01T01_AC01_E13W02_R01_Z1_A1.mat');
[~,idx] = ismember({'Pitch','Surge'},labels);

basin_time = channels(:,1);
basin_pitch = channels(:,idx(1));
basin_surge = channels(:,idx(2));

% Load OpenFAST results
load('C:\Umaine Google Sync\GitHub\SHARK\OpenFAST\Simulations\FOCAL_C4\Elsevier_LC02\FOCAL_C4_FAST_Results.mat');

fast_time = sim_results.Time;
fast_pitch = sim_results.PtfmPitch;
fast_surge = sim_results.PtfmSurge;

%% ----- Plot Comparisons ----- %%
fig = figure
ax = gca; hold on; box on;
ax = setAxesInfo(ax);
xlim([1000,1400])
plot(basin_time,basin_pitch*(pi/180),'DisplayName','Experiment','LineWidth',1.5,'Color','Black')
plot(fast_time,fast_pitch*(pi/180),'DisplayName','OpenFAST','LineWidth',1.5)

xlabel('Time [s]')
ylabel('Angle [rads]')
title({'Tuned OpenFAST Model vs','Experimental Pitch Angle'})
legend

exportgraphics(fig,'Figures/FOCAL_Tuned_FAST_Pitch.pdf','ContentType','vector')

fig = figure;
ax = gca; hold on; box on;
ax = setAxesInfo(ax);
xlim([1000,1400])

plot(basin_time,basin_surge,'DisplayName','Experiment','LineWidth',1.5,'Color','Black')
plot(fast_time,fast_surge,'DisplayName','OpenFAST','LineWidth',1.5)

xlabel('Time [s]')
ylabel('Displacement [m]')
title({'Tuned OpenFAST Model vs','Experimental Surge Displacement'}')
legend

exportgraphics(fig,'Figures/FOCAL_Tuned_FAST_Surge.pdf','ContentType','vector')