function xc = EKF(dt, x, RHSFunc, JacoFunc, MeasFunc, Q, R, y, varargin)

% xp = EKF(dt, x, RHSFunc, JacoFunc, Q, R, y) performs a state update using an
% extended Kalman filter (EKF) algorithm. State update is first performed
% with a runge-kutta 4th order integration scheme on the update function
% RHSFunc. Next, the linear matrices are computed using the supplied
% JacoFunc and correction is done using the measurements supplied by
% MeasFunc and corresponding system measurements (y). 
% 
% xp = EKF(dt, x, RHSFunc, JacoFunc, MeasFunc, Q, R, y) does the integration including an input
% F to RHSFunc.
%
% Inputs:
% ---------
% % dt - time step [s]
%
% x - current state vector (n x 1)
%
% RHSFunc - handle to function for computing state derivative given the
% current state. Can be linear or non-linear but must take the form
% RHSFunc(x, u) where x is the state and u is the input.
%
% JacoFunc - handle to function which compute the system jacobians from the
% current state. Must take the form [A, H] = JacoFunc(dt, x, F) where A is the
% linear system dynamics matrix and H is the linear output matrix
%
% MeasFunc - handle to function which computes model outputs from state
% vector and input, of the form MeasFunc(x, F)
%
% Q - starting state noise covariance matrix
%
% R - starting measurement noise covariance matrix
%
% F (optional) - input terms for the current time step formatted accoring
% to RHSFunc requirements
%
% Written By: Ian Ammerman
% Last Modified: 8/29/24

%% ----- Parse Inputs ----- %%
p = inputParser;

% Required inputs
addRequired(p, 'dt')
addRequired(p, 'x')
addRequired(p, 'RHSFunc', @validateFunctionHandle)
addRequired(p, 'MeasFunc', @validateFunctionHandle)
addRequired(p, 'JacoFunc', @validateFunctionHandle)
addRequired(p, 'Q', @ismatrix)
addRequired(p, 'R', @ismatrix)
addRequired(p, 'y')

% Optional arguments
addOptional(p, 'F', 0, @validateForcingHandle)

% Parse
parse(p, dt, x, RHSFunc, MeasFunc, JacoFunc, Q, R, y, varargin{:})

% Distribute inputs
dt = p.Results.dt;
x = p.Results.x;
RHSFunc = p.Results.RHSFunc;
JacoFunc = p.Results.JacoFunc;
MeasFunc = p.Results.MeasFunc;
Qinit = p.Results.Q;
Rinit = p.Results.R;
y = p.Results.y;
F = p.Results.F;

%% ----- Persistent Variables ----- %%
persistent Qc Rc P S init_flag

if isempty(init_flag)
    Qc = Qinit;
    Rc = Rinit;
    P = zeros(size(Qinit));
    init_flag = 1;
end

%% ----- Perform Prediction Step ----- %%
xp = RK4(dt, x, RHSFunc, F);

%% ----- Get Linear System ----- %%
[A, H] = JacoFunc(dt, xp, F);

%% ----- Get Model Measurements ----- %%
yp = MeasFunc(xp, F);

%% ----- Perform Update Step ----- %%
% Update P
P = A*P*A' + Qc;

% Compute measurement residual
res = y - yp;

% Compute Kalman gain
S = H*P*H' + Rc;
K = P*H'/S;

% Update state
xc = xp + K*res;

% Update P
P = (eye(size(K*H)) - K*H)*P + 1e-6 * eye(size(P));

end

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