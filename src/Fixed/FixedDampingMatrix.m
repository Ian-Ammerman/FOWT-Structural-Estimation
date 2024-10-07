function C = FixedDampingMatrix(K)

persistent btable ttable init_flag

if isempty(init_flag)
    [btable,ttable] = load5MWDistributed;
    init_flag = 1;
end

% Define structural damping ratios
zt = 0.01;
zb = 0.17;

%% ----- Tower structural damping ----- %%
% omegat = 1.1889 * 2 * pi; % tower natural frequency in rads/s
omegat = 0.32 * 2 * pi; % tower natural frequency in rads/s
mt = trapz(ttable.Ht, ttable.TMassDen .* ttable.mode1.^2);
ct = 2*zt*mt*omegat;
% ct = 2*(zt/omegat)*K(1,1);

%% ----- Tower Aerodynamic Damping ----- %%
% Fd = trapz(0.5*rho*D*phi^2(h)q^2Cd(h)
% Fd = 0.5*rho*q^2 * trapz(Ht, D*phi^2(h))
% ct_aero = 0.5 * 1.225 * trapz(ttable.Ht, ttable.TwrDiam .* ttable.mode1.^3);
ct_aero = 0;
%% ----- Blade structural damping ----- %%
omegab = 0.75 * 2*pi; % blade flapwise natural frequency in rads/s
cb1 = ((2*zb)/(omegab)) * K(2,2);
cb2 = ((2*zb)/(omegab)) * K(3,3);
cb3 = ((2*zb)/(omegab)) * K(4,4);

%% ----- Drivetrain Damping ----- %%
P0 = 5*10^6; % rated mechanical power
Q0 = 12.1*(pi/30);
% cdt = -(P0/Q0^2);
cdt = 6215000;

% Form damping matrix
C = diag([ct + ct_aero, cb1, cb2, cb3]);

end