function xdot = FixedTrinityRHS(x, F)

persistent btable ttable init_flag

if isempty(init_flag)
    [btable,ttable] = load5MWDistributed;

    init_flag = 1;
end

pos = x(1:length(x)/2);
vel = x(1+length(x)/2:end);

azimuth = F.azimuth + [0, 2*pi/3, 4*pi/3];
omega = F.omega;

M = FixedMassMatrix(btable, ttable);
K = FixedStiffnessMatrix(btable, ttable, azimuth, omega);
C = FixedDampingMatrix(K);

% State accelerations
a = M \ (F.forces - C*vel - K*pos);

% Compile xdot
xdot = [vel;a];

