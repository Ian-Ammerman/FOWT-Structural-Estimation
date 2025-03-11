clear all; close all; clc;

load('Data/Elsevier_LC02/5MW_OC4Semi_WSt_WavesWN_FAST_Results.mat')

%% ----- Compare Nacelle Pitch Velocity to True Pitch Velocity ----- %%
fig = figure;
% axis padded
ax = gca; hold on; box on;
ax = setAxesInfo(ax);
xlim([220,280])

plot(sim_results.Time,sim_results.QD_P,'DisplayName','True','LineWidth',1.5,'Color','Black')
plot(sim_results.Time,sim_results.NcIMURVys*(pi/180),'DisplayName','IMU','LineWidth',1)

xlabel('Time [s]')
ylabel('Angular Velocity [rads/s]')
title('Platform Pitch Velocity Comparison')
legend

exportgraphics(fig,'IMU_Pitch_Validation.pdf','ContentType','vector')