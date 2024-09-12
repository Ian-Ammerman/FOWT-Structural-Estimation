function [xplus, kbar] = RK4(dt, x, RHSFunc, varargin)

% xp = RK4(dt, x, RHSFunc) performs 4th order Runge-Kutta integration on
% the function RHSFunc given the current state x and time step dt, assuming
% free vibration.
%
% xp = RK4(dt, x, RHSFunc, F) performs the integration including the
% external force/input term F, which must match the required input for
% RHSFunc.
%
% Inputs:
% ----------
% dt - time step [s]
%
% x - current state vector (n x 1)
%
% RHSFunc - handle to function for computing state derivative given the
% current state. Can be linear or non-linear but must take the form
% RHSFunc(x, u) where x is the state and u is the input.
%
% F (optional) - input terms for the current time step formatted accoring
% to RHSFunc requirements
%
% Written By: Ian Ammerman
% Last Modified: 8/29/24

%% ----- Parse Inputs ----- %%
p = inputParser;

% Required inputs (time step, state)
addRequired(p, 'dt')
addRequired(p, 'x', @isvector)
addRequired(p, 'RHSFunc', @validateFunctionHandle)

% Optional input (system inputs/forces)
addOptional(p, 'F', struct(), @validateForcingHandle)

% Parse inputs
parse(p, dt, x, RHSFunc, varargin{:})

% Store inputs
dt = p.Results.dt;
x = p.Results.x;
RHSFunc = p.Results.RHSFunc;
F = p.Results.F;

%% ----- Do State Update (Runge-Kutta 4th Order Integration) ----- %%
% Compute factors
k1 = RHSFunc(x, F);
k2 = RHSFunc(x + dt*(k1/2), F);
k3 = RHSFunc(x + dt*(k2/2), F);
k4 = RHSFunc(x + dt*k3, F);

% Get mean acceleration
kbar = (1/6)*(k1 + 2*k2 + 2*k3 + k4);

% Update state
xplus = x + dt*kbar;

%% ----- Helper Functions ----- %%
function validateFunctionHandle(input)
    if ~isa(input, 'function_handle')
        error('Input must be a function handle.');
    end
end

function validateForcingHandle(input)
    if ~isstruct(input) && ~isvector(input) && ~isscalar(input)
        error('Input must be a structure, vector, or scalar.');
    end
end

end