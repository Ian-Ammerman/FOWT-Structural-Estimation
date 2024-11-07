function y = NREL5MW_Measurement(x,F)

Rotor = F.rotor;
Tower = F.tower;

% Initialize output vector
y = zeros(2,1);

% Platform pitch angle
y(1,1) = x(1);

% Tower velocity
% y(2,1) = x(8);

% Tower base bending moment
y(2,1) = Tower.TwFAStif(1) * Tower.ddmode1(1) * x(2) * 10^-8;

% Root bending load
y(3,1) = pi*Rotor.FlpStff(1) * Rotor.BladeddFmode1(1) * x(3) * 10^-6;
y(4,1) = pi*Rotor.FlpStff(1) * Rotor.BladeddFmode1(1) * x(4) * 10^-6;
y(5,1) = pi*Rotor.FlpStff(1) * Rotor.BladeddFmode1(1) * x(5) * 10^-6;