function xdot = RHS(x, F, fM, fC, fK)

% Separate position and velocity states
pos = x(1:length(x)/2);
vel = x(length(x)/2+1:end);

% Form mass, stiffness, damping matrices
M = fM(F.azimuth,F.tower,F.rotor);
K = fK(F.azimuth,F.omega,F.tower,F.rotor);
C = fC(K,M);

% Form RHS forces
rhs = F.forces - C*vel - K*pos;

% Solve for accelerations
acc = M\rhs;

% Output state derivative
xdot = [vel;acc];