%% Code for solving Min Enegy consumption for a altitude difference for Vertical Take-off

%% Use of ode45 commands + fminbnd for calculatin optimum.
%% THIS FUNCTION DOES NOT RETURN DISTRIBUTIONS, BUT ONLY MORE PRECISE OPTIMALS
clear all
close all
%% Set deltaH (altitude difference)
global deltaH
deltaH = 100; %[m]
global D
D     = 0.8124; % Propeller Diameter [m]
global N_eng;
N_eng = 2; % Number of engines
global rho
rho   = 1.2133; % Air density [kg/m^3]

RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;
%% Optimisation calculation
% RPS that give minimum Energy consumption

n_opt = fminbnd(@(n) accel_n_cte(n), 55, rps_max)


%% Definition of the ode45 that solves the whole problem

function [E] = accel_n_cte (n)

%% Define time interpolation
Tlin = 100; % Steps for time
tspan = linspace(0,30,Tlin);      % Interpolate time of integration
y0 = [0 0];                       % Initial conditions. H(0) = 0; V(0) = 0
TOL = 1e-10;                      % Define tolerance

%% Definition of options
optionsolver = odeset('AbsTol',TOL,'RelTol',TOL,'Events',@eventsH);

[t,y] = ode45(@(t,y)fun(t,y,n),tspan,y0,optionsolver);

global D
global N_eng
global rho
%% Electric
tau = 0.2;
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency    
%% Power model

CP0 = 0.0261518307541734;
CP1 = 0.0473735972985378;
CP2 = -0.16267474946046;
CP3 = 0.0247028469343899;

    V_vec = interp1(linspace(0,1,numel(y(:,2))),y(:,2),linspace(0,1,Tlin)); % For speed
    t_vec = interp1(linspace(0,1,numel(t)),t,linspace(0,1,Tlin));           % For time
    
    
    %% Energy & Power Calculations
    %% Define CP as a vector and not as a function to increase calculation speed within the loop
    P_vec =  (CP3*V_vec.^3./(n.^3*D^3)+CP2*V_vec.^2./(n.^2*D^2)+CP1*V_vec./(n.*D)+CP0).*rho*N_eng*n.^3*D^5;
    
    E     = trapz(t_vec,P_vec)/((1-tau)*eta_m);

end

%% Definition of function to be solved

function dy = fun(t,y,n)
% Some parameters are not required within this function, but they are
% included anyway. Will delete in a further version of the code.

global D
global N_eng
global rho

%% Aerodynamic model
alpha_V = -90*pi/180; % [º] Because of Vertical TO


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


% Full model w.r.t angle of attack

CD = @(alpha) p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + p_CD_ac(4)*alpha.^3 ...
    + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);

%% Propulsive Model

CT0 =  0.0735880531010883;
CT1 = -0.0311758018412727;
CT2 = -0.249744726429543; 




%% Loads and weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
%W     = mTOW*g;

S = 0.5008; % Reference surface [m^2]




%% Thrust
CT = @(V,n) CT2*V^2/(n^2*D^2)+CT1*V/(n*D)+CT0;

T = @(V,n) CT(V,n)*N_eng*rho*n^2*D^4;




%% Define system of equations

dy = zeros(2,1);
dy(1) = y(2); % dh/dt = V
dy(2) = T(y(2),n)/mTOW - g - 1/2*rho*y(2).^2*S*CD(alpha_V)/mTOW; % m*dV/dt = T-W-D

end


function [position,isterminal,direction] = eventsH(t,y)
global deltaH
position = y(1)-deltaH; % The value that we want to be zero
isterminal = 1;           % Halt integration
direction = 0;            % The zero can be approached from either direction
end


