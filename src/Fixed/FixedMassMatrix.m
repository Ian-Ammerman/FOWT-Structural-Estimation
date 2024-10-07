function M = FixedMassMatrix(btable, ttable)

%% ----- Flapwise Only ----- %%
persistent mtb mtt mbb mNH IDT Irotor IgenLSS init_flag



if isempty(init_flag)

    % Nacelle mass
    mNH = 110000 + 240000;

    % Tower mass term
    mtt = 3*trapz(btable.Rb, btable.BMassDen) + trapz(ttable.Ht, ttable.TMassDen .* ttable.mode1.^2) + mNH;

    % Blade mass term
    mbb = trapz(btable.Rb, btable.BMassDen .* btable.Fmode1.^2);

    % Tower-blade mass coupling terms
    mtb = trapz(btable.Rb, btable.BMassDen .* btable.Fmode1);

    % Rotor and drivetrain
    N = 97;
    Irotor = 35447886 + 115926;
    Igen = 534.116;
    IgenLSS = Igen * N^2;
    IDT = Irotor + IgenLSS;

    init_flag = 1;
end

% Form mass matrix
% M = [mtt, mtb, mtb, mtb,                0;
%      mtb, mbb,   0,   0,                0;
%      mtb,   0, mbb,   0,                0;
%      mtb,   0,   0, mbb,                0;
%        0,   0,   0,   0, Irotor + IgenLSS];

M = [mtt, mtb, mtb, mtb;
     mtb, mbb,   0,   0;
     mtb,   0, mbb,   0;
     mtb,   0,   0, mbb];













