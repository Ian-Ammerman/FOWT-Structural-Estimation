 % Basic OpenFAST driver for new turbine testing models
clear all; close all; clc;

%% ---------- USER INPUTS ---------- %%
home_dir = 'C:\Umaine Google Sync\GitHub\SHARK';

% Simulation folder for outputs
sim_folder = 'Turbulent_13mps_Monhegan_NoWave';

% OpenFAST Model to Run
model = '5MW_Baseline_Fixed';

% OpenFAST version (name of corresponding folder in 'bin')
version = 'v3_5_3';
bin_name = 'openfast_x64.exe';

%% ---------- END USER INPUTS ---------- %%

%% ---------- RUN SIMULATION ---------- %%
% Ensure starting from home dir
if ~strcmp(pwd,home_dir)
    cd(home_dir)
end

% Directory Definitions 
bin_dir = sprintf('OpenFAST/bin/%s/%s',version,bin_name);
model_dir = sprintf('OpenFAST/%s',model);
fst_path = sprintf('%s/%s.fst', model_dir, model);
sim_dir = sprintf('OpenFAST/Simulations/%s/%s',model,sim_folder);

% Go home
cd(home_dir)

% Check for simulation directory
if ~exist(sim_dir,'dir')
    answer = questdlg(sprintf('Simulation directory "%s" not found. Would you like to create it?', sim_folder),...
                              'Simulation Folder Creation', ...
                              'Yeah, sure', 'OMG NO!', 'Yes');
    switch answer
        case 'Yeah, sure'
            mkdir(sim_dir)
        case 'OMG NO!'
            return
    end
end

% Run OpenFAST Simulation
name = sprintf('"%s/%s" "%s/%s"', home_dir, bin_dir, home_dir, fst_path);
[status,results] = system(name,'-echo');

%% ---------- RELOCATE OUTPUT FILES ---------- %%
disp('---------- Relocating Output Files -----------')
cd(sim_dir)

filetypes = {'.out','.ech','.sum','.outb','.yaml'};
for i = 1:length(filetypes)
    try
        movefile(sprintf('../../../%s/*%s',model,filetypes{i}));
        fprintf('Output files of type %s relocated. \n',filetypes{i});
    catch
        fprintf('No files of type %s detected. \n',filetypes{i});
    end
end

%% ----- Process Outputs into Matlab Format ----- %%
% Process output files & format as structure (assumes .out format)
[sim_results,units] = readFastTabular(sprintf('%s.out', model));

% Save results structure named after model used
save_name = sprintf('%s_FAST_Results.mat',model);
save(save_name,'sim_results');

% Relocate
fclose('all');

%% ----- Plot Some Results ----- %%
figure
gca; hold on; box on;
plot(sim_results.Time,sim_results.RotSpeed,'DisplayName','OpenFAST');
title('Rotor Speed [RPM]')
xlabel('Time [s]')

figure
gca; hold on; box on;
plot(sim_results.Time, sim_results.TTDspFA)
title('Tower-Top Displacement [m]')
xlabel('Time [s]')

figure
gca; hold on; box on;
plot(sim_results.Time, sim_results.BldPitch1)
title('Blade Pitch Command [deg]')
xlabel('Time [s]')