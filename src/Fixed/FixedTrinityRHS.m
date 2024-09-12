function xdot = FixedTrinityRHS(x, F)

global btable ttable

pos = x(1:length(x)/2);
vel = x(1+length(x)/2:end);

azimuth = F.azimuth + [0, 2*pi/3, 4*pi/3];
omega = F.omega;

M = FixedMassMatrix(btable, ttable);
K = FixedStiffnessMatrix(btable, ttable, azimuth, omega);
C = FixedDampingMatrix(K);

% State accelerations
a = M \ (F.forces - C*vel - K*pos);

% Jerk (constant acceleration)
j = 0;

% Compile xdot
% xdot = [vel;a;j];
xdot = [vel;a];

