function K = RotorStiffnessMatrix(btable, azimuth, omega)

persistent init_flag kg_con ke ktt kgr_con

% Initialize constant values
if isempty(init_flag)

    % Elastic stiffness
    ke = trapz(btable.Rb, btable.FlpStff .* btable.ddFmode1.^2);

    % Pre-compute gravity term for centrifugal and gravitational stiffnesses
    mu1 = zeros(length(btable.Rb),1);
    for i = 1:length(btable.Rb)-1
        mu1(i) = trapz(btable.Rb(i:end), btable.BMassDen(i:end).*btable.Rb(i:end));
    end
    mu1(end) = mu1(end-1);

    mu2 = zeros(length(btable.Rb),1);
    for i = 1:length(btable.Rb)-1
        mu2(i) = trapz(btable.Rb(i:end), btable.BMassDen(i:end));
    end
    mu2(end) = mu2(end-1);

    % Centrifugal stiffening constant term
    kg_con = trapz(btable.Rb, mu1 .* btable.dFmode1.^2);
    
    % Gravitational stiffening constant term
    kgr_con = -0.5 * 9.806 * trapz(btable.Rb, mu2 .* btable.dFmode1.^2);

    init_flag = 1;
end

% Centrifugal stiffness
kg = omega^2 * kg_con;

% Gravitational stiffening
kgr = cos(azimuth) * kgr_con;

% Combine for final terms
kb = kgr + ke + kg;

%% ----- Form Stiffness Matrix ----- %%
K = diag([kb(1), kb(2), kb(3)]);
end