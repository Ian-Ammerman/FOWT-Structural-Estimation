function K = NREL5MW_StiffnessMatrix(azimuth, omega, tower, rotor, floater)

persistent init_flag kg_con ke kdt kpp ktt kgr_con
N = 97;

% azimuth = azimuth + [0, 2*pi/3, 4*pi/3];

% Initialize constant values
if isempty(init_flag)

    % Elastic stiffness
    ke = trapz(rotor.BladeStations, rotor.FlpStff .* rotor.BladeddFmode1.^2);

    % Pre-compute gravity term for centrifugal and gravitational stiffnesses
    mu1 = zeros(length(rotor.BladeStations),1);
    for i = 1:length(rotor.BladeStations)-1
        mu1(i, 1) = trapz(rotor.BladeStations(i:end), rotor.BMassDen(i:end).*rotor.BladeStations(i:end));
    end
    mu1(end) = mu1(end-1);

    mu2 = zeros(length(rotor.BladeStations),1);
    for i = 1:length(rotor.BladeStations)-1
        mu2(i, 1) = trapz(rotor.BladeStations(i:end), rotor.BMassDen(i:end));
    end
    mu2(end) = mu2(end-1);
    
    % Centrifugal stiffening constant term
    kg_con = trapz(rotor.BladeStations, mu1 .* rotor.BladeddFmode1.^2);
    
    % Gravitational stiffening constant term
    kgr_con = -0.5 * 9.806 * trapz(rotor.BladeStations, mu2 .* rotor.BladeddFmode1.^2);
    
    % Tower stiffness
    ktt = trapz(tower.Ht, tower.TwFAStif .* tower.ddmode1.^2);
    
    % Platform stiffness
    kpp = 1.75*(floater.p_HST + floater.p_mooring_stiffness); % Hydrostatic + mooring stiffness

    init_flag = 1;
end

% Centrifugal stiffness
kg = omega^2 * kg_con;

% Gravitational stiffening
kgr = cos(azimuth) * kgr_con;

% Combine for final terms
kb = kgr + ke + kg;

%% ----- Form Stiffness Matrix ----- %%
K = diag([kpp, ktt, kb(1), kb(2), kb(3)]);
end