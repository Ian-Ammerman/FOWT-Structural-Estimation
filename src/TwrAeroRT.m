function F = TwrAeroRT(u,qT,alpha)

persistent ttable z rho init_flag

if isempty(init_flag)
    [~,ttable] = load5MWDistributed;
    z = ttable.Ht + 10;
    rho = 1.225;
    init_flag = 1;
end

%% ----- Project Wind w/ Power Law ----- %%
uinf = powerLaw(u,90,z,alpha) - qT*flip(ttable.mode1);

%% ----- Compute Drag Force ----- %%
F = 0.5*rho*trapz(ttable.Ht,ttable.TwrDiam.*uinf.^2.*ttable.mode1);

%% ----- Helper Function ----- %%
function u = powerLaw(uref,zref,z,alpha)
    u = uref.*(z./zref).^alpha;
end

end