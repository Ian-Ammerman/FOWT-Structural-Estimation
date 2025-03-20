% function TurbSimDriver(input_case)

input_case = 'H90m_13mps_Monhegan';

%% ----- Define Directories ----- %%
% Fixed
home_dir = 'C:\Users\ianam\Documents\GitHub\FOWT-Structural-Estimation\OpenFAST\Turbsim';
exe_rel = 'bin\TurbSim_x64.exe';
exe_global = sprintf('%s/%s', home_dir, exe_rel);

% Dynamic (from input)
inp_file = sprintf('%s/%s.inp',input_case, input_case);
inp_global = sprintf('%s/%s',home_dir, inp_file);

%% ----- Execute TurbSim ----- %%
[status, cmdout] = system(sprintf('"%s" "%s"', exe_global, inp_global), '-echo');