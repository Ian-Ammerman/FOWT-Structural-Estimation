function [J, H, B] = FixedTrinityJacobian(dt, x, F)

persistent btable ttable sigma_con init_flag

if isempty(init_flag)
    [btable,ttable] = load5MWDistributed;
    % Constant for bending stress calculation
    sigma_con = (3.542/2) * trapz(btable.Rb, btable.Rb .* btable.Fmode1 .* btable.FlpStff);
    init_flag = 1;
end

%% ----- Prepare Operating Point ----- %%
% Compute OP
xop = x;

%% ----- Perform State Linearization ----- %%
% Define perturbation
dx = 0.0001*[1,-1];
for i = 1:length(x)
    for j = 1:2
        % Perturb x
        x(i) = xop(i) + dx(j);

        % Compute new x
        xtemp(:,j) = RK4(dt, x, @FixedTrinityRHS, F);

        % Reset x
        x = xop;        
    end

    % Compute jacobian
    Jx(:,i) = (xtemp(:,2) - xtemp(:,1))/(dx(2) - dx(1));
end

% Output linearization
J = Jx;

%% ----- Perform Output Linearization ----- %%
% % Define perturbation
% dx = 0.0001*[1,-1];
% for i = 1:length(x)
%     for j = 1:2
%         % Perturb x
%         x(i) = xop(i) + dx(j);
% 
%         % Compute new x
%         ytemp(:,j) = FixedTrinityMeasurement(x,F);
% 
%         % Reset x
%         x = xop;        
%     end
% 
%     % Compute jacobian
%     Jy(:,i) = (ytemp(:,2) - ytemp(:,1))/(dx(2) - dx(1));
% end

H = [0,         0,         0,         0, 1, 0, 0, 0;
     0, sigma_con,         0,         0, 0, 0, 0, 0;
     0,         0, sigma_con,         0, 0, 0, 0, 0;
     0,         0,         0, sigma_con, 0, 0, 0, 0];

%% ----- Perform Input Linearization ----- %%
B = 0;



end