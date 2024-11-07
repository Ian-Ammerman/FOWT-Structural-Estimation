function Rotor = loadNREL5MWRotor

%% ----- Define Blade ----- %%
% Define constant parameters
Lb = 63;

% Load in distributed parameters
blade_aero = readmatrix('OpenFAST\5MW_Baseline\AeroData\NREL_5MW_Blade.csv');
blade = readtable("OpenFAST\5MW_Baseline\StructData\blade_distributed.txt");

% Define radial stations (minus hub)
blade.Rb = Lb * blade.BlFract;

% Define blade fundamental mode polynomials
F1_coeffs = [-2.2555, 4.7131, -3.2452, 1.7254, 0.0622, 0, 0];
E1_coeffs = [-0.6952, 2.376, -3.5772, 2.5337, 0.3627, 0, 0];

% Adjust coefficients to be function of x alone
idx = [6, 5, 4, 3, 2, 1, 0];
F1_coeffs = F1_coeffs./Lb.^idx;
E1_coeffs = E1_coeffs./Lb.^idx;

dF1_coeffs = polyder(F1_coeffs);
dE1_coeffs = polyder(E1_coeffs);

ddF1_coeffs = polyder(dF1_coeffs);
ddE1_coeffs = polyder(dE1_coeffs);

% Evaluate mode shapes and derivatives
blade.Fmode1 = polyval(F1_coeffs, blade.Rb);
blade.Emode1 = polyval(E1_coeffs, blade.Rb);

blade.dFmode1 = polyval(dF1_coeffs, blade.Rb);
blade.dEmode1 = polyval(dE1_coeffs, blade.Rb);

blade.ddFmode1 = polyval(ddF1_coeffs, blade.Rb);
blade.ddEmode1 = polyval(ddE1_coeffs, blade.Rb);

blade.EIo = blade.FlpStff .* cosd(blade.StrcTwst) + blade.EdgStff .* sind(blade.StrcTwst);

%% ----- Define Rotor Properties ----- %%
% Rotor properties
Rotor.RH = 1.5;                         % hub radius [m]
Rotor.R = max(blade.Rb);                % blade radius [m]
Rotor.B = 3;                            % number of blades
Rotor.CLTable = readmatrix('OpenFAST\5MW_Baseline\AeroData\NREL_5MW_Airfoil_CL_Lookup.csv');
Rotor.CDTable = readmatrix('OpenFAST\5MW_Baseline\AeroData\NREL_5MW_Airfoil_CD_Lookup.csv');
Rotor.AlphaTable = readmatrix('OpenFAST\5MW_Baseline\AeroData\NREL_5MW_Airfoil_Alpha_Lookup.csv');

% Blade properties
Rotor.blade = blade_aero;               % blade distributed aerodynamic properties
Rotor.BladeStations = blade.Rb;         % radial station locations [m]
Rotor.BMassDen = blade.BMassDen;        % blade linear mass density
Rotor.BladeFmode1 = blade.Fmode1;       % 1st flapwise bending mode shape normalized at tip
Rotor.BladedFmode1 = blade.dFmode1;     % derivative of 1st flapwise mode shape
Rotor.BladeddFmode1 = blade.ddFmode1;   % second derivative of 1st flapwise mode shape

Rotor.FlpStff = blade.FlpStff;

% Dynamic properties
Rotor.RPM = [0];
Rotor.dqb = [0,0,0];