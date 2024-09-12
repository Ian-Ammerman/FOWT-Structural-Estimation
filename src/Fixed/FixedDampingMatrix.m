function C = FixedDampingMatrix(K)

% Define structural damping ratios
zt = 0.025;
% zb = 0.005;
zb = 0.17;

% Tower structural damping
omegat = 0.32 * 2*pi; % tower natural frequency in rads/s
ct = ((2*zt)/(omegat)) * K(1,1);

% Blade structural damping
omegab = 0.75 * 2*pi; % blade flapwise natural frequency in rads/s
cb1 = ((2*zb)/(omegab)) * K(2,2);
cb2 = ((2*zb)/(omegab)) * K(3,3);
cb3 = ((2*zb)/(omegab)) * K(4,4);

% Form damping matrix
C = diag([ct, cb1, cb2, cb3]);

end