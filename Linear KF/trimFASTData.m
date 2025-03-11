function t = trimFASTData(dt)

% Load in data
load('FOCAL_C4_FAST_Results.mat');
to_load = {'Time','PtfmPitch','NcIMUTVxs','NcIMURVys','TTDspFA','RotThrust'};
data = zeros(length(sim_results.Time),length(to_load));
for i = 1:length(to_load)
    data(:,i) = sim_results.(to_load{i});
end

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
t = array2table(newdata,'VariableNames',to_load);

% Get velocity vectors
chans = {'PtfmPitch'};
newchans = {'PtfmPitchVel'};
dt = mean(diff(t.Time));
for i = 1:length(chans)
    vals = gradient(t.(chans{i}),dt);
    
    newbar = mean(vals);
    vals = lowpass(vals-newbar,1,1/dt,'ImpulseResponse','iir');
    t.(newchans{i}) = vals + newbar;
end

%% ----- Trim Time ----- %%
tmin = 800;
tmax = 1200;

idxmin = dsearchn(t.Time,tmin);
idxmax = dsearchn(t.Time,tmax);

newt = table;
for i = 1:width(t)
    varName = t.Properties.VariableNames{i};
    newt.(varName) = t.(varName)(idxmin:idxmax);
end

t = newt;
t.Time = t.Time - min(t.Time);