function C = RotorDampingMatrix(K)

% Define structural damping ratios
zb = 0.005;
zb = 0.17;

% Blade structural damping
omegab = 0.68 * 2*pi; % blade flapwise natural frequency in rads/s
cb1 = ((2*zb)/(omegab)) * K(1,1);
cb2 = ((2*zb)/(omegab)) * K(2,2);
cb3 = ((2*zb)/(omegab)) * K(3,3);

% Form damping matrix
C = diag([cb1, cb2, cb3]);

end