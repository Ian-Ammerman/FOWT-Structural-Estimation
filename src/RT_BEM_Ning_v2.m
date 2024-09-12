function forces = RT_BEM_Ning_v2(wind_vector,dqb,RPM,blade_pitch)

persistent blade CL_vals CD_vals btable Alpha_vals R RH B rho r beta dr c at Fmode1 ne sp init_flag

if isempty(init_flag)
    %Load blade node data
    %Columns: 1- radial position (m), 2 - twist (deg), 3 - dr (m)
    %3 - chord (m), %4 - %airfoil # (-)
    blade = readmatrix('C:\Umaine Google Sync\Masters Working Folder\Multibody FOWT Model\Aero Data\NREL_5MW_Blade.csv');
    
    % Load CL/CD lookup information
    CL_vals = readmatrix('C:\Umaine Google Sync\Masters Working Folder\Multibody FOWT Model\Aero Data\NREL_5MW_Airfoil_CL_Lookup.csv');
    CD_vals = readmatrix('C:\Umaine Google Sync\Masters Working Folder\Multibody FOWT Model\Aero Data\NREL_5MW_Airfoil_CD_Lookup.csv');
    Alpha_vals = readmatrix('C:\Umaine Google Sync\Masters Working Folder\Multibody FOWT Model\Aero Data\NREL_5MW_Airfoil_Alpha_Lookup.csv');

    % Load in structural information
    [btable, ~] = load5MWDistributed;

    %Set inputs (These inputs are for the NREL 5 MW reference wind turbine)
    R = 63; %Blade radius (m)
    RH = 1.5; %Hub radius (m)
    B = 3;
    
    rho = 1.225; %Density of air (kg/m^3)
    
    %Separate out blade information
    r = blade(:,1); %Node radial position
    beta = (blade(:,2))*pi/180; %Node twist (converted to radians)
    dr = blade(:,3); %Element width dr
    c = blade(:,4); %Node chord
    at = blade(:,5); %Airfoil #

    % Interpolate mode shape to match aerodynamics table
    Fmode1 = pchip(btable.Rb, btable.Fmode1, r);
    
    %Determine number of elements
    ne = length(r);
    
    %Compute local solidity
    sp = (1/(2*pi))*3.*c./r;

    init_flag = 0;
end

a = 0;
ap = 0;

% Angular speed in rads/s
omega = RPM * (pi/30);

[thrust,torque,rt,rn,flap,edge,rf,re] = deal([0;0;0]);
for j = 1:3

    % Local inflow wind speed
    Uinf = wind_vector(:,j);

    %Initialize total thrust and total torque
    Ttotal = 0; %Thrust (N)
    Qtotal = 0; %Torque (N-m)
    Ltotal = 0;
    Dtotal = 0;
    FlapTotal = 0;
    
    % Initialize incremental thrust, torque, lift, and drag
    [dTr,dQr,dLr,dDr] = deal(zeros(size(Uinf)));

 
    % Loop over blade elements
    for i = ne:-1:1   
        
        %% ----- Estimate Inflow Velocities ----- %%
        Vx = Uinf(i) - dqb(j)*Fmode1(i);
        Vy = omega*r(i);

        %% ----- Get Force Coefficients from Sectional Theory ----- %%
        if abs(Vy) >= 0
            %Compute inflow angle (rad)
            phi = atan(Vx/Vy);
    
            %Estimate angle of attack in degrees
            alpha = (180/pi)*(phi-beta(i)) - (180/pi)*blade_pitch(j);
    
            % Obtain lift and drag coefficients
            CL = interp1(Alpha_vals,CL_vals(:,at(i)),alpha);
            CD = interp1(Alpha_vals,CD_vals(:,at(i)),alpha);
            
            % Compute normal force coefficient
            cx = CL*cos(phi) + CD*sin(phi);
    
            % Compute tangential force coefficient
            cy = CL*sin(phi) - CD*cos(phi);
            
            %% ----- Apply Prandtl's Correction for Hub/Tip Loss Factors ----- %%
            % Compute loss factors
            ftip = (B/2)*((R - r(i))/(r(i)*abs(sin(phi))));
            Ftip = (2/pi)*acos(exp(-ftip));
    
            fhub = (B/2)*((r(i)-RH)/(RH*abs(sin(phi))));
            Fhub = (2/pi)*acos(exp(-fhub));
    
            % Combine loss factors
            F = Ftip*Fhub;
    
            % Define convinience parameters
            k = (sp(i)*cx)/(4*F*sin(phi)^2);
            kp = (sp(i)*cy)/(4*F*sin(phi)*cos(phi));
    
            %% ----- Determine Solution Region & Induction Factors ----- %%
            % Axial induction factor
            if phi > 0 && k <= (2/3) % Momentum region
                a = k/(1+k);
            elseif phi > 0 && k > (2/3) % Empirical region
                g1 = 2*F*k - ((10/9) - F);
                g2 = 2*F*k - F*((4/3)-F);
                g3 = 2*F*k - ((25/9)-2*F);
                
                if g3 == 0 
                    a = 1 - 1/(2*sqrt(g2));
                else
                    a = (g1 - sqrt(g2))/g3;
                end
            elseif phi < 0 && k > 1 % Propeller brake region
                a = k/(k-1);
            elseif phi < 0 && k <= 0 % No possible solution here
                a = 0;
            end
            
            % Tangential induction factor
            ap = kp/(1-kp);
    
            %% ----- Define Residual Function ----- %%
            Rphi = @(p) (sin(p)/(1-a)) - (Vx/Vy)*(cos(phi)/(1+ap));
                
            %% ----- Solve Residual Function for New Phi ----- %%
            if Rphi(eps)*Rphi(pi/2) < 0
                % fplot(Rphi,[eps,pi/2]); drawnow
                phi_new = fzero(Rphi,[eps,pi/2]);
            elseif Rphi(-pi/4)*Rphi(-eps) < 0
                % fplot(Rphi,[-pi/4,-eps]); drawnow
                phi_new = fzero(Rphi,[-pi/4,-eps]);
            else
                % fplot(Rphi,[pi/2,pi-eps]); drawnow
                phi_new = fzero(Rphi,[pi/2,pi-eps]);
            end
    
            % Update phi
            phi = phi_new;
        else
            phi = pi/2;
        end

        %% ----- Determine New Force Coefficients from Sectional Theory ----- %%
        %Estimate angle of attack in degrees
        alpha = (180/pi)*(phi-beta(i)) - (180/pi)*blade_pitch(j);

        % Obtain lift and drag coefficients
        CL = interp1(Alpha_vals,CL_vals(:,at(i)),alpha);
        CD = interp1(Alpha_vals,CD_vals(:,at(i)),alpha);

        % Recompute inflow velocities
        Vx = Uinf(i)*(1-a);
        Vy = RPM*(2*pi/60)*r(i)*(1+ap);

        %Compute Ustar
        Ustar = sqrt(Vx^2 + Vy^2); 
        
        %% ----- Compute Aerodynamic Loads ----- %%
        % Lift force
        dL = 0.5*rho*(Ustar^2)*CL*c(i)*dr(i);
        Ltotal = Ltotal + dL;
        dLr(i) = dL*r(i);

        % Drag force
        dD = 0.5*rho*(Ustar^2)*CD*c(i)*dr(i);
        Dtotal = Dtotal + dD;
        dDr(i) = dD*r(i);

        % Rotor thrust
        dT = dL*cos(phi) + dD*sin(phi);
        Ttotal = Ttotal+dT; 
        dTr(i) = dT*r(i);

        % Rotor torque
        dQ = (dL*sin(phi)-dD*cos(phi))*r(i);
        Qtotal = Qtotal+dQ;
        dQr(i) = dQ*r(i);

        % Flapwise force
        dFlap = (dL*cosd(alpha) + dD*sind(alpha))*Fmode1(i);
        FlapTotal = FlapTotal + dFlap;
    
    end
    thrust(j,1) = Ttotal;
    torque(j,1) = Qtotal;
    rt(j,1) = sum(dQr)/Qtotal;
    rn(j,1) = sum(dTr)/Ttotal;
    flap(j,1) = FlapTotal;
    edge(j,1) = Dtotal;
    rf(j,1) = sum(dLr)/Ltotal;
    re(j,1) = sum(dDr)/Dtotal;
end

check_flag = 0;

if check_flag == 1
    forces = thrust;
elseif check_flag == 0
    forces.Thrust = thrust;
    % forces.Torque = torque;
    forces.FlapwiseForce = flap;
    % forces.EdgewiseForce = edge;
    % forces.Rt = rt;
    % forces.Rn = rn;
    % forces.Rf = rf;
    % forces.Re = re;
    % forces.Fn = thrust;
    % forces.Ft = torque./rt;
end