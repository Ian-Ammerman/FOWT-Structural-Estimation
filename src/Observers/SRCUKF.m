function xc = SRCUKF(dt, x, RHSFunc, MeasFunc, Q, R, y, varargin)

%% ----- Parse Inputs ----- %%
p = inputParser;

% Required inputs
addRequired(p, 'dt')
addRequired(p, 'x')
addRequired(p, 'RHSFunc', @validateFunctionHandle)
addRequired(p, 'MeasFunc', @validateFunctionHandle)
addRequired(p, 'Q', @ismatrix)
addRequired(p, 'R', @ismatrix)
addRequired(p, 'y')

% Optional arguments
addOptional(p, 'F', 0, @isstruct)

% Parameters
addParameter(p, 'Weight', 1)

% Parse
parse(p, dt, x, RHSFunc, MeasFunc, Q, R, y, varargin{:})

% Distribute inputs
dt = p.Results.dt;
x = p.Results.x;
RHSFunc = p.Results.RHSFunc;
MeasFunc = p.Results.MeasFunc;
Qinit = p.Results.Q;
Rinit = p.Results.R;
y = p.Results.y;
F = p.Results.F;
sigma_weight = p.Results.Weight;

%% ----- Persistent Variables ----- %%
persistent Qc Rc L W1 W2 eta Sk Xp Sqk Srk init_flag

if isempty(init_flag)
    Qc = Qinit;
    Rc = Rinit;

    L = length(x); % length of state vector

    % SCRKF Tuning Parameters
    eta = sigma_weight;

    % Matrix square roots
    Sqk = chol(Qc);
    Srk = chol(Rc);
    Sk = chol(Qc);

    % Init sigma point matrices
    Xp = zeros(L, 2*L);

    W1 = 1/(2*L);
    W2 = sqrt(W1);
    
    init_flag = 1;
end

%% ----- Generate Sigma Points ----- %%
X = repmat(x, [1,2*L]) + eta*[Sk, -Sk];

%% ----- Predict State ----- %%
% Update sigma points
for i = 1:size(X,2)
    Xp(:,i) = RK4(dt, X(:,i), RHSFunc, F);
end

% Re-combine for mean prediction
xp = W1 * sum(Xp,2);

% Update cholesky factor in time
[Sk_minus, ~] = qr([W2 * (Xp - repmat(xp, [1,2*L])), Sqk]);

%% ----- Generate More Sigma Points ----- %%
X2 = repmat(xp, [1,2*L]) + eta*[Sk_minus, -Sk_minus];

%% ----- Compute Outputs ----- %%
% Get measurements
for i = 1:size(X2,2)
    Yp(:,i) = MeasFunc(X2(:,i), F);
end

% Mean measurement
yp = W1 * sum(Yp, 2);

% Update cholesky factor in time
[Syk, ~] = qr([W2 * (Yp - repmat(yp, [1,2*L])), Srk]);

%% ----- Compute Cross-Covariance & Do State Update ----- %%
% Compute cross-covariance
Pxy = W1 * (X2 - repmat(xp, [1,2*L])) * (Yp - repmat(yp, [1,2*L]))';

% Compute Kalman gain
K = (Pxy/Syk')/Syk;

% Compute residual
res = y - yp;

% Correct state vector
xc = xp + K*res;

%% ----- Posterior Covariance Update ----- %%
% Define convinience parameter
Zk = W2 * ((X2 - repmat(xp, [1,2*L])) - K*(Yp - repmat(yp, [1,2*L])));

% Update cholesky factor
[Sk, ~] = qr([Zk, K*Srk]);


%% ----- Helper Functions ----- %%
function validateFunctionHandle(input)
    if ~isa(input, 'function_handle')
        error('Input must be a function handle.');
    end
end

end

















