% clear all; close all; clc;
% dt = 0.01;
function t = trimData(dt)

% LC01
% load('Data/FOCAL_4.2_M02D01T03_AC01_E02W01_R01_Z1_A1.mat');

% LC02
load('Data/FOCAL_4.2_M02D01T01_AC01_E13W02_R01_Z1_A1.mat');

% LC03
% load('Data/FOCAL_4.2_M02D01T05_AC01_E23W03_R01_Z1_A1.mat');

to_load = [1,5:10,12,13,24:29,30,40];

tags = labels(to_load);
data = channels(:,to_load);

% Update time vector to match dt
time = data(:,1);
newtime = [min(time):dt:max(time)]';
newdata(:,1) = newtime;

for i = 2:size(data,2)
    bar = mean(data(:,i));
    data(:,i) = lowpass(data(:,i)-bar,1,1/mean(diff(time)),'ImpulseResponse','iir');
    data(:,i) = data(:,i) + bar;
    newdata(:,i) = pchip(time,data(:,i),newtime);
    % newdata(:,i) = smoothdata(newdata(:,i),'sgolay');
end

% Form table
t = array2table(newdata,'VariableNames',tags);

% Get velocity vectors
chans = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
newchans = {'SurgeVel','SwayVel','HeaveVel','RollVel','PitchVel','YawVel'};
dt = mean(diff(t.Time));
for i = 1:length(chans)
    vals = gradient(t.(chans{i}),dt);
    
    newbar = mean(vals);
    vals = lowpass(vals-newbar,1,1/dt,'ImpulseResponse','iir');
    t.(newchans{i}) = vals + newbar;
end

%% ----- Trim Time ----- %%
tmin = 700;
tmax = 1300;

idxmin = dsearchn(t.Time,tmin);
idxmax = dsearchn(t.Time,tmax);

newt = table;
for i = 1:width(t)
    varName = t.Properties.VariableNames{i};
    newt.(varName) = t.(varName)(idxmin:idxmax);
end

t = newt;
t.Time = t.Time - min(t.Time);

% Add tower deflection
EI = 2499900000000.00;
dv2dx = 0.000147858426311794;
t.TTDspFA = t.towerBotMy./(EI*dv2dx);