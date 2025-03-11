clear all; close all; clc;

%% ----- Load in LC02 Results ----- %%
% OpenFAST baseline
load('Data/Elsevier_LC02/5MW_OC4Semi_WSt_WavesWN_FAST_Results.mat');

% Estimator results
load('Data/EKF_LC02.mat');
load('Data/EKF_wFBG_LC02.mat');
load('Data/UKF_LC02.mat');
load('Data/UKF_wFBG_LC02.mat');


%% ----- Plot Comparison ----- %%
% Platform pitch velocity
trange = [150,180];

fig = figure;
ax = gca; hold on; box on;
ax = setAxesInfo(ax);
xlim(trange);

plot(sim_results.Time,sim_results.QD_P,'DisplayName','OpenFAST','LineWidth',1.5,'Color','Black')
plot(EKF(:,1),EKF(:,7),'DisplayName','EKF','LineWidth',1.5,'LineStyle','--','Color','Red');
plot(EKF_wFBG(:,1),EKF_wFBG(:,7),'DisplayName','EKF w/FBG','LineWidth',1.5,'LineStyle',':','Color','Red');
% plot(UKF(:,1),UKF(:,7),'DisplayName','UKF','LineWidth',1.5,'LineStyle','-.','Color','Blue');
% plot(UKF_wFBG(:,1),UKF_wFBG(:,7),'DisplayName','UKF w/FBG','LineWidth',1.5,'LineStyle',':','Color','Blue');

xlabel('Time [s]')
ylabel({'Platform Pitch Angular','Velocity [rads/s]'})
% title('Platform Pitch Velocity')
legend

% exportgraphics(fig,'Figures/LC02_EKF_UKF_Pitch_Velocity.pdf','ContentType','vector')

% Tower-top displacement
trange = [150,200];

fig = figure;
ax = gca; hold on; box on;
ax = setAxesInfo(ax);
xlim(trange);

plot(sim_results.Time,sim_results.TTDspFA,'DisplayName','OpenFAST','LineWidth',1.5,'Color','Black')
plot(EKF(:,1),EKF(:,3),'DisplayName','EKF','LineWidth',1.5,'LineStyle','--','Color','Red');
plot(EKF_wFBG(:,1),EKF_wFBG(:,3),'DisplayName','EKF w/FBG','LineWidth',1.5,'LineStyle',':','Color','Red');
% plot(UKF(:,1),UKF(:,3),'DisplayName','UKF','LineWidth',1.5,'LineStyle','-.','Color','Blue');
% plot(UKF_wFBG(:,1),UKF_wFBG(:,3),'DisplayName','UKF w/FBG','LineWidth',1.5,'LineStyle',':','Color','Blue');

xlabel('Time [s]')
ylabel('Tower 1^{st} Mode Displacement [m]')
% title('Tower 1^{st} Mode Displacement')
legend

% exportgraphics(fig,'Figures/LC02_EKF_UKF_Tower_Comparison.pdf','ContentType','vector')

% Blade 1 flapwise displacement
trange = [160,180];

fig = figure;
ax = gca; hold on; box on;
ax = setAxesInfo(ax);
xlim(trange);

plot(sim_results.Time,sim_results.Q_B1F1,'DisplayName','OpenFAST','LineWidth',1.5,'Color','Black')
plot(EKF(:,1),EKF(:,4),'DisplayName','EKF','LineWidth',1.5,'LineStyle','--','Color','Red');
plot(EKF_wFBG(:,1),EKF_wFBG(:,4),'DisplayName','EKF w/FBG','LineWidth',1.5,'LineStyle',':','Color','Red');
% plot(UKF(:,1),UKF(:,4),'DisplayName','UKF','LineWidth',1.5,'LineStyle','-.','Color','Blue');
% plot(UKF_wFBG(:,1),UKF_wFBG(:,4),'DisplayName','UKF w/FBG','LineWidth',1.5,'LineStyle',':','Color','Blue');

xlabel('Time [s]')
ylabel('Blade 1 1^{st} Mode Displacement [m]')
% title('Blade 1 1^{st} Mode Displacement')
legend

% exportgraphics(fig,'Figures/LC02_EKF_UKF_Blade_Displacement_Comparison.pdf','ContentType','vector')


% Blade 1 flapwise velocity
trange = [160,180];

fig = figure;
ax = gca; hold on; box on;
ax = setAxesInfo(ax);
xlim(trange);

plot(sim_results.Time,sim_results.QD_B1F1,'DisplayName','OpenFAST','LineWidth',1.5,'Color','Black')
% plot(EKF(:,1),EKF(:,9),'DisplayName','EKF','LineWidth',1.5,'LineStyle','--','Color','Red');
% plot(EKF_wFBG(:,1),EKF_wFBG(:,9),'DisplayName','EKF w/FBG','LineWidth',1.5,'LineStyle',':','Color','Red');
plot(UKF(:,1),UKF(:,9),'DisplayName','UKF','LineWidth',1.5,'LineStyle','-.','Color','Blue');
plot(UKF_wFBG(:,1),UKF_wFBG(:,9),'DisplayName','UKF w/FBG','LineWidth',1.5,'LineStyle',':','Color','Blue');

xlabel('Time [s]')
ylabel('Blade 1 1^{st} Modal Velocity [m/s]')
% title('Blade 1 1^{st} Modal Velocity')
legend

% exportgraphics(fig,'Figures/LC02_EKF_UKF_Blade_Velocity_Comparison.pdf','ContentType','vector')