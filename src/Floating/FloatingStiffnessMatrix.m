function K = FloatingStiffnessMatrix(azimuth, omega)

persistent init_flag kg_con ke kdt kpp kpt_con ktt kgr_con
N = 97;

azimuth = azimuth + [0, 2*pi/3, 4*pi/3];

% Initialize constant values
if isempty(init_flag)

    [btable,ttable] = load5MWOC4SemiDistributed;

    % Elastic stiffness
    ke = trapz(btable.Rb, btable.FlpStff .* btable.ddFmode1.^2);

    % Pre-compute gravity term for centrifugal and gravitational stiffnesses
    mu1 = zeros(length(btable.Rb),1);
    for i = 1:length(btable.Rb)-1
        mu1(i, 1) = trapz(btable.Rb(i:end), btable.BMassDen(i:end).*btable.Rb(i:end));
    end
    mu1(end) = mu1(end-1);

    mu2 = zeros(length(btable.Rb),1);
    for i = 1:length(btable.Rb)-1
        mu2(i, 1) = trapz(btable.Rb(i:end), btable.BMassDen(i:end));
    end
    mu2(end) = mu2(end-1);
    
    % Centrifugal stiffening constant term
    kg_con = trapz(btable.Rb, mu1 .* btable.dFmode1.^2);
    
    % Gravitational stiffening constant term
    kgr_con = -0.5 * 9.806 * trapz(btable.Rb, mu2 .* btable.dFmode1.^2);
    
    % Tower stiffness
    ktt = trapz(ttable.Ht, ttable.TwFAStif .* ttable.ddmode1.^2);
    
    % Platform stiffness
    kpp = 1.75*(3.8032E8 + 8.6E7); % Hydrostatic + mooring stiffness
    
    % Drivetrain stiffness
    Irotor = 35447886 + 115926;
    Igen = 534.116;
    IDT = Irotor + Igen * N^2;
    
    kdt = 0.2^2 * IDT;

    kpt_con = trapz(btable.Rb, btable.Rb.*btable.BMassDen);

    init_flag = 1;
end

% Centrifugal stiffness
kg = omega^2 * kg_con;

% Gravitational stiffening
kgr = cos(azimuth) * kgr_con;

% Combine for final terms
kb = kgr + ke + kg;

%% ----- Form Stiffness Matrix ----- %%
K = diag([kpp, ktt, kb(1), kb(2), kb(3), kdt]);
end