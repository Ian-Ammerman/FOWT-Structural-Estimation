function C = NREL5MW_DampingMatrix(K, M)

% Define structural damping ratios
zp = 0.01;
zt = 0.01;
zb = 0.17;

% Pitch linear damping
cp = 2 * zp * sqrt(K(1,1) * M(1,1));

% Tower structural damping
omegat = 1.889 * 2*pi; % tower natural frequency in rads/s
% ct = 2*zt*omegat*M(2,2);
ct = ((2*zt)/(omegat)) * K(2,2);

% Blade structural damping
omegab = 0.75 * 2*pi; % blade flapwise natural frequency in rads/s
cb1 = ((2*zb)/(omegab)) * K(3,3);
cb2 = ((2*zb)/(omegab)) * K(4,4);
cb3 = ((2*zb)/(omegab)) * K(5,5);

% Drivetrain damping
P0 = 5*10^6; % rated mechanical power
Q0 = 12.1*(pi/30);
cdt = -(P0/Q0^2);

% Form damping matrix
C = diag([cp, ct, cb1, cb2, cb3, cdt]);

end