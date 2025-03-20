function xc = SRCUKF(dt, x, RHSFunc, MeasFunc, Q, R, y, F, eta)

%% ----- Persistent Variables ----- %%
persistent Qc Rc L W1 W2 Sk Xp Sqk Srk init_flag

if isempty(init_flag)
    Qc = Q;
    Rc = R;

    L = length(x); % length of state vector

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
Yp = zeros(size(Rc,1),2*L);
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

end

















