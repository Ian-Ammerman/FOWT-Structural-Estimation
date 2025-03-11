function Tower = loadNREL5MWTower

%% ----- Define Tower ----- %%
% Define constant parameters
Ht = 87.6;


% Read in distributed parameters
tower = readtable("C:\Umaine Google Sync\Masters Working Folder\1 - OpenFAST\Models/5MW_Baseline/StructData/OC4Semi_Tower_Distributed.txt");
tower_aero = readtable("C:\Umaine Google Sync\Masters Working Folder\1 - OpenFAST\Models/5MW_Baseline/AeroData/tower_aero.txt");

% Define height stations
tower.Ht = tower.HtFract * Ht;

% Interpolate aerodynamic values along tower height
tower.TwrDiam = pchip(tower_aero.TwrElev,tower_aero.TwrDiam,tower.Ht);
tower.Cd = ones(size(tower.Ht));

% Define fundamental mode
idx = [6, 5, 4, 3, 2, 1, 0];
T1_coeffs = [0.4082, -1.5035, 1.8042, -0.8622, 1.1533, 0, 0];
T1_coeffs = T1_coeffs./Ht.^idx;

dT1_coeffs = polyder(T1_coeffs);

ddT1_coeffs = polyder(dT1_coeffs);

% Evaluate mode shape and derivatives
tower.mode1 = polyval(T1_coeffs, tower.Ht);
tower.dmode1 = polyval(dT1_coeffs, tower.Ht);
tower.ddmode1 = polyval(ddT1_coeffs, tower.Ht);

Tower = table2struct(tower,"ToScalar",true);
Tower.NacHubMass = 110000 + 240000;
Tower.HH = 87.6;