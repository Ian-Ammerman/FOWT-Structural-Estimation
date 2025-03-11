function xdot = RK4RHS(t,x)

%% ----- Define Dynamics Model ----- %%
% Mass matrix
m11 = 48410870305.1658;
m22 = 1278084.80047176;
m12 = 185459490.208370;
m21 = 185459490.208370;
M = [m11, m12;
     m21, m22];

% Damping matrix
c1 = 197713958.925553;
c2 = 20701.7415674758;
C = [c1,  0;
      0, c2];

% Stiffness matrix
k1 = 4.0848e+09;  % Hydrostatic + mooring stiffness
k2 = 3.8235e+06;  % Tower modal stiffness
K = [k1,  0;
      0, k2];

%% ----- Get States ----- %%
% size(x)
pos = x(1:length(x)/2,1);
vel = x(length(x)/2 + 1:end,1);

%% ----- Solve for Acceleration ----- %%
rhs = -C*vel - K*pos;
acc = M\rhs;

xdot = [vel;acc];


