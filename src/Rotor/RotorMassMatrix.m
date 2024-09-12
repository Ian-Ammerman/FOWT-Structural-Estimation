function M = RotorMassMatrix(btable)

%% ----- Flapwise Only ----- %%
persistent mbb init_flag

if isempty(init_flag)

    % Blade mass term
    mbb = trapz(btable.Rb, btable.BMassDen .* btable.Fmode1.^2);

    init_flag = 1;
end

% Form mass matrix
M = diag([mbb, mbb, mbb]);