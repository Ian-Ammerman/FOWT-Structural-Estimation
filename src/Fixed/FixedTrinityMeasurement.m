function y = FixedTrinityMeasurement(x, F)

persistent btable ttable kappa_norm init_flag

if isempty(init_flag)
    [btable,ttable] = load5MWDistributed;

    init_flag = 1;
end

% Tower-top velocity
y(1,1) = [0, 0, 0, 0, 1, 0, 0, 0]*x;

% Root bending load
y(2,1) = pi*btable.FlpStff(1) * btable.ddFmode1(1) * x(2) * 10^-6;
y(3,1) = pi*btable.FlpStff(1) * btable.ddFmode1(1) * x(3) * 10^-6;
y(4,1) = pi*btable.FlpStff(1) * btable.ddFmode1(1) * x(4) * 10^-6;