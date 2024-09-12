function xdot = RotorRHS(x, F)

global btable

pos = x(1:3);
vel = x(4:6);

azimuth = F.RotAz + [0, 2*pi/3, 4*pi/3];
omega = F.omega;

M = RotorMassMatrix(btable);
K = RotorStiffnessMatrix(btable, azimuth, omega);
C = RotorDampingMatrix(K);

a = M \ (F.forces - C*vel - K*pos);
xdot = [vel;a];