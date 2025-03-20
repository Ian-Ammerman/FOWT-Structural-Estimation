function xc = KF(dt, x, u, meas, Q, R, sys)

% Written by: Ian Ammerman
% Last Modified: 11/27/2024
% 
% xc = KF(dt,x,u,y,Q,R,sys) simulates a linear Kalman filter applied to the
% system object sys using a discrete time formulation.
%
% Inputs:
% ----------
% dt - time step [s]
% x  - current state
% u  - current input (applied forces)
% meas - vector of measurement signals
% Q  - state covariance matrix
% R  - measurement covariance matrix
% sys - state-space system object

persistent F B C D P Ts dsys init_flag

%% ----- Initialize Persistent Variables ----- %%
if isempty(init_flag)
    % Check for Discrete System 
    if ~isdt(sys)
        dsys = c2d(sys,dt);
    else
        dsys = sys;
    end
    
    % Extract System Information
    [F,B,C,D,Ts] = ssdata(dsys);

    % Initialize Filter
    P = Q;

    init_flag = 0;
end

%% ----- Peform Prediction Step ----- %%
[x,P] = predict(x,P,F,Q,B,u);

%% ----- Perform Update Step ----- %%
[xc,P,~] = update(C,P,R,meas,x);

%% ----- Helper Functions ----- %%
% Prediction (Labbe, 2020, pg 212)
function [x,P] = predict(x,P,F,Q,B,u)
    x = F*x + B*u; %predict states
    P = F*P*F' + Q; %predict process covariance
end

% Update
function [x,P,K] = update(H,P,R,z,x)
    S = H*P*H' + R; % Project system uncertainty into measurement space & add measurement uncertainty
    K = P*H'*inv(S);
    y = z-H*x; % Error term
    x = x+K*y;
    KH = K*H;
    P = (eye(size(KH))-KH)*P;
end

end