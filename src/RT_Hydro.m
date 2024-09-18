function F = RT_Hydro(dt, x, xdot)

persistent init_flag Kt Ainf vel_store

if isempty(init_flag)
    %% ----- Prepare for Convolution Calculation ----- %%
    % Radiation memory effect
    load('C:\Umaine Google Sync\Masters Working Folder\Multibody FOWT Model\FAST_Model\HydroData\Ainf.mat');
    load('C:\Umaine Google Sync\Masters Working Folder\Multibody FOWT Model\FAST_Model\HydroData\Aomega.mat');
    load('C:\Umaine Google Sync\Masters Working Folder\Multibody FOWT Model\FAST_Model\HydroData\Bomega.mat');
    load('C:\Umaine Google Sync\Masters Working Folder\Multibody FOWT Model\FAST_Model\HydroData\Omega.mat');
    
    % Get periods and frequencies
    wRs = w;
    wHz = wRs/(2*pi);
    periods = 1./wHz;
    
    % Define the sampling frequency (f_s)
    f_s = 2 * max(wHz); 
    
    % Number of points (length of the IFFT data)
    N = length(w);
    
    % Sampling interval
    T_s = 1 / f_s;
    
    % Time vector
    time_vector = (0:N-1) * T_s;
    
    % Compute K(t)
    for i = 1:length(time_vector)
        cos_term = reshape(cos(wRs*time_vector(i)),[1,1,length(time_vector)]);
        Kt(:,:,i) = (2/pi)*trapz(w,Bw.*cos_term,3);
    end
    
    % Interpolate to new time vector
    newtime = [0:0.0125:max(time_vector)];
    Ktnew = zeros(size(Kt,1),size(Kt,2),length(newtime));
    for i = 1:size(Kt,1)
        for j = 1:size(Kt,2)
            vals = squeeze(Kt(i,j,:));
            newvals = pchip(time_vector,vals,newtime);
            Ktnew(i,j,:) = reshape(newvals,[1,1,length(newtime)]);
        end
    end
    
    % Only take last 60 seconds
    idx = dsearchn(newtime',60);
    Kt = Ktnew(:,:,1:idx);
    N = size(Kt,3);

    vel_store = zeros(6,N);

    init_flag = 1;
end

%% ----- Load in New Velocity/Acceleration Values ----- %%
% Assign current acceleration/velocity values
vtemp = [0, 0, 0, 0, x(1), 0]';
acc = xdot(1);

% Update stored values
vel_store = [vtemp, vel_store(:,1:end-1)];

%% ----- Compute Radiation Loads ----- %%
f_rad = zeros(6,1);
for j = 1:size(vel_store,2)
    f_rad = f_rad + Kt(:,:,j)*vel_store(:,j);
end
f_rad = -Ainf*acc - dt*f_rad;

%% ----- Output Hydrodynamic Forces ----- %%
F.PtfmPitch = f_rad(5);