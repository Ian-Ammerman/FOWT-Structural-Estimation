function C = FloatingDampingMatrix(K, M, azimuth, omega)

persistent init_flag btable ttable xmux

azimuth = azimuth + [0, 2*pi/3, 4*pi/3];

if isempty(init_flag)
    [btable, ttable] = load5MWDistributed;

    xmux = trapz(btable.Rb, btable.Rb.*btable.BMassDen);

    init_flag = 1;
end

% Define structural damping ratios
zp = 0.01;
zt = 0.01;
zb = 0.17;

% Pitch linear damping
cp = 2 * zp * sqrt(K(1,1) * M(1,1));

% Tower structural damping
omegat = 0.32 * 2*pi; % tower natural frequency in rads/s
ct = ((2*zt)/(omegat)) * K(2,2);

ctplus = sum(omega * sin(azimuth) * xmux);

% Blade structural damping
omegab = 0.68 * 2*pi; % blade flapwise natural frequency in rads/s
cb1 = ((2*zb)/(omegab)) * K(3,3);
cb2 = ((2*zb)/(omegab)) * K(4,4);
cb3 = ((2*zb)/(omegab)) * K(5,5);

% Drivetrain damping
P0 = 5*10^6; % rated mechanical power
Q0 = 12.1*(pi/30);
cdt = -(P0/Q0^2);

% Form damping matrix
C = diag([cp, ct, cb1, cb2, cb3, cdt]);

C = [    cp, ctplus,   0,   0,    0,   0;
     ctplus,     ct,   0,   0,    0,   0;
          0,      0, cb1,   0,    0,   0;
          0,      0,   0, cb2,    0,   0;
          0,      0,   0,   0,  cb3,   0;
          0,      0,   0,   0,    0, cdt];

end