function y = FloatingTrinityMeasurement(x,F)

persistent btable ttable kappa_norm init_flag

if isempty(init_flag)
    [btable,ttable] = load5MWOC4SemiDistributed;

    init_flag = 1;
end

% Initialize output vector
y = zeros(2,1);

% Platform pitch velocity
y(1,1) = x(7);

% Tower velocity
y(2,1) = x(8);
% 
% % Root bending load
% y(3,1) = pi*btable.FlpStff(1) * btable.ddFmode1(1) * x(3) * 10^-6;
% y(4,1) = pi*btable.FlpStff(1) * btable.ddFmode1(1) * x(4) * 10^-6;
% y(5,1) = pi*btable.FlpStff(1) * btable.ddFmode1(1) * x(5) * 10^-6;