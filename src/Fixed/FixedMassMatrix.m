function M = FixedMassMatrix(btable, ttable)

%% ----- Flapwise Only ----- %%
persistent mtb mtt mbb mNH init_flag

if isempty(init_flag)

    % Nacelle mass
    mNH = 110000 + 240000;

    % Tower mass term
    mtt = 3*trapz(btable.Rb, btable.BMassDen) + trapz(ttable.Ht, ttable.TMassDen .* ttable.mode1.^2) + mNH;

    % Blade mass term
    mbb = trapz(btable.Rb, btable.BMassDen .* btable.Fmode1.^2);

    % Tower-blade mass coupling terms
    mtb = trapz(btable.Rb, btable.BMassDen .* btable.Fmode1);

    init_flag = 1;
end

% Form mass matrix
M = [mtt, mtb, mtb, mtb;
     mtb, mbb,   0,   0;
     mtb,   0, mbb,   0;
     mtb,   0,   0, mbb];