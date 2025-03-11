clear all; close all; clc;

%% ----- Load in FOCAL Data ----- %%
load('FOCAL_4.2_M02D01T01_AC01_E13W02_R01_Z1_A1.mat')
dt = 0.01;
t = trimData(dt);

g = 9.806;      % gravitational acceleration
rt = 145.23;    % height from SWL to nacelle accelerometer
rb = 15.25;     % height from SWL to tower base accelerometer

%% ----- Get Pitch Acceleration ----- %%
t.SurgeAcc = gradient(t.SurgeVel,dt);
t.HeaveAcc = gradient(t.HeaveVel,dt);
t.PitchAcc = gradient(t.PitchVel*(pi/180),dt);

ttvel = cumtrapz(t.Time,highpass(t.accelNacelleAx + 9.806*t.Pitch*(pi/180),0.05,1/dt));
ttvel_adj = ttvel - rt*t.PitchVel*(pi/180);

figure
gca; hold on
plot(t.Time,t.accelNacelleAx + 9.806*t.Pitch*(pi/180),'DisplayName','Raw')
% plot(t.Time,highpass(t.accelNacelleAx + 9.806*t.Pitch*(pi/180),0.05,1/dt),'DisplayName','HP Filtered')
% plot(t.Time,lowpass(t.accelNacelleAx + 9.806*t.Pitch*(pi/180),0.05,1/dt,'ImpulseResponse','iir'),'DisplayName','LP Filtered')
legend

% % Plot dominant x accelerations at nacelle
% figure
% gca; hold on; box on;
% plot(t.Time,t.SurgeAcc,'DisplayName','Surge');
% plot(t.Time,t.HeaveAcc,'DisplayName','Heave');
% plot(t.Time,rt*t.PitchAcc,'DisplayName','Pitch Tangential');
% legend
% 
% % Plot sum
% figure
% gca; hold on; box on;
% plot(t.Time,bandpass(rt*t.PitchAcc,[0.05,0.15],1/dt),'DisplayName','Sum')
% plot(t.Time,bandpass(t.accelNacelleAx + 9.806*t.Pitch*(pi/180),[0.05,0.15],1/dt),'DisplayName','Nacelle Measured')
% % plot(t.Time,highpass(t.towerBotAccX + 9.806*t.Pitch*(pi/180),0.05,1/dt),'DisplayName','Platform Measured')
% legend

figure
gca; hold on; box on;
plot(t.Time,ttvel,'DisplayName','Integrated')
plot(t.Time,ttvel_adj,'DisplayName','Isolated')
plot(t.Time,highpass(ttvel,0.15,1/dt),'DisplayName','Filtered')
% plot(t.Time,rt*t.PitchVel*(pi/180),'DisplayName','Measured')
legend

% %% ----- Compare to Nacelle Acceleration ----- %%
% ar = rt*t.PitchAcc;
% 
% figure
% gca; hold on; box on;
% plot(t.Time,ar-9.806*t.Pitch*(pi/180),'DisplayName','from Pitch (linear)');
% plot(t.Time,ar-9.806*sind(t.Pitch),'DisplayName','from Pitch (nonlinear)')
% % plot(t.Time,t.accelNacelleAx,'DisplayName','Measured')
% legend

% %% ----- Compare Nacelle Acceleration to Platform Acceleration ----- %%
% figure
% gca; hold on; box on;
% plot(t.Time,t.accelNacelleAx,'DisplayName','Nacelle')
% plot(t.Time,t.towerBotAccX,'DisplayName','Tower Base')
% plot(t.Time,t.accelNacelleAx - ar,'DisplayName','Adjusted Nacelle')
% legend

% %%
% figure
% gca; hold on;
% plot(t.Time,t.SurgeVel)
% 
% figure
% gca; hold on; 
% plot(t.Time,t.Surge,'DisplayName','Processed')
% plot(channels(:,1),channels(:,5),'DisplayName','Original')
% % plot(t.Time,lowpass(t.Surge,1.5,1/0.01,'ImpulseResponse','iir'),'DisplayName','Lowpass')
% legend
% 
% figure
% gca; hold on;
% plot(channels(:,1),channels(:,24),'DisplayName','Original')
% plot(t.Time,t.accelNacelleAx,'DisplayName','Processed','LineWidth',1)

% %% ----- Plot FFT of Acceleration ----- %%
% af = myFFT(t.accelNacelleAx,1/0.01);
% afbar = myFFT(rMean(t.accelNacelleAx),1/0.01);
% afbarhp = myFFT(highpass(rMean(t.accelNacelleAx),0.05,1/0.01,'ImpulseResponse','iir')+mean(t.accelNacelleAx),1/0.01);
% 
% figure
% gca; hold on; box on;
% plot(af(:,1),af(:,2),'DisplayName','Default')
% plot(afbar(:,1),afbar(:,2),'DisplayName','De-meaned')
% plot(afbarhp(:,1),afbarhp(:,2),'DisplayName','Highpass')
% xlim([0,0.2])
% ylim([0,0.1])
% 
% %%
% 
% ah = gradient(t.HeaveVel,0.01);
% 
% bar = mean(t.accelNacelleAx);
% hp = highpass(t.accelNacelleAx-bar, 0.05, 1/0.01,'ImpulseResponse','iir') + bar;
% 
% 
% 
% figure
% gca; hold on; 
% plot(t.Time,-9.806*t.Pitch*(pi/180),'DisplayName','From Acceleration')
% plot(t.Time,t.accelNacelleAx,'DisplayName','Processed')
% plot(t.Time,highpass(rMean(t.accelNacelleAx),0.08,1/0.01,'ImpulseResponse','iir')+mean(t.accelNacelleAx),'DisplayName','Highpass')
% % plot(channels(:,1),channels(:,24),'DisplayName','Original')
% legend