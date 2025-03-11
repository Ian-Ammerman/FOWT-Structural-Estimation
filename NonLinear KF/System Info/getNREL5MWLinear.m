function [A,C] = getNREL5MWLinear(dt,x,u,azimuth,omega)

persistent init_flag rotor tower floater

if isempty(init_flag)

    % Load system information
    rotor = loadNREL5MWRotor;
    tower = loadNREL5MWTower;
    floater = loadNREL5MWFloater;

    init_flag = 1;
end

%% ----- Define Dynamics Matrices ----- %%
% Load mass matrix
M = NREL5MW_MassMatrix(azimuth,tower,rotor,floater);

% Stiffness matrix
K = NREL5MW_StiffnessMatrix(azimuth,omega,tower,rotor,floater);

% Damping matrix
C = NREL5MW_DampingMatrix(K,M);

% Form state-space system
A = [zeros(size(M)),eye(size(M));
     -M\K, -M\C];

B = [zeros(size(M));-M\eye(size(M))];

%% ----- Define Output Function ----- %%
% No blade measurements
C = [1,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,1,0,0,0,0;
     0,0,0,0,0,0,1,0,0,0];

% % Blade IMUs
% C = [1,0,0,0,0,0,0,0,0,0;
%      0,0,0,0,0,1,0,0,0,0;
%      0,0,0,0,0,0,1,0,0,0;
%      0,0,0,0,0,0,0,1,0,0;
%      0,0,0,0,0,0,0,0,1,0;
%      0,0,0,0,0,0,0,0,0,1];

% % Blade FBG
% C = [1,0,0,0,0,0,0,0,0,0;
%      0,0,0,0,0,1,0,0,0,0;
%      0,0,0,0,0,0,1,0,0,0;
%      0,0,1,0,0,0,0,0,0,0;
%      0,0,0,1,0,0,0,0,0,0;
%      0,0,0,0,1,0,0,0,0,0];

D = zeros(size(C,1),size(B,2));

%% ----- Discretize ----- %%
sys = ss(A,B,C,D);
dsys = c2d(sys,dt);

[A,B,C,D] = ssdata(dsys);

end

