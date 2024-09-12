function y = FixedTrinityMeasurement(x, F)

persistent btable ttable kappa_norm init_flag

if isempty(init_flag)
    [btable,ttable] = load5MWDistributed;

    % Constant for bending stress calculation
    % kappa_norm = trapz(btable.Rb, btable.FlpStff .* btable.ddFmode1);
    kappa_norm = btable.FlpStff(1) * btable.ddFmode1(1);

    init_flag = 1;
end

% Compute mass, damping, and stiffness matrix
M = FixedMassMatrix(btable,ttable);
K = FixedStiffnessMatrix(btable,ttable,F.azimuth + [0, 2*pi/3, 4*pi/3],F.omega);
C = FixedDampingMatrix(K);

% Tower-top velocity
y(1,1) = [0, 0, 0, 0, 1, 0, 0, 0]*x;

% Root bending load
% y(2,1) = K(2,2) * x(2);

% % Tower-top acceleration
% [~, acc_full] = RK4(0.0125, x, @FixedTrinityRHS, F);
% y(2,1) = acc_full(1);

% % Blade flapwise bending load
% sigma_flap = kappa_norm * [x(2), x(3), x(4)]';
% y(2:4,1) = sigma_flap;