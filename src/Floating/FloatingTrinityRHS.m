function xdot = FloatingTrinityRHS(x, F)

% Separate position and velocity states
pos = x(1:6);
vel = x(7:12);

% Form mass, stiffness, damping matrices
M = FloatingMassMatrix(F.azimuth,F.omega);
K = FloatingStiffnessMatrix(F.azimuth,F.omega);
C = FloatingDampingMatrix(K,M);

% Form RHS forces
rhs = F.forces - C*vel - K*pos;

% Solve for accelerations
acc = M\rhs;

% Output state derivative
xdot = [vel;acc];