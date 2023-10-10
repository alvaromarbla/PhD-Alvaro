clc
clear all
close all

%% Calculate cruise from V1 to V2 given a certain RPS. The angle of attack and speed will be time dependant.
%  This problem will be solved in an numerical manner using Differential Algebraic Equations Solver.
%  The solution requires the integration of a Differential Algebraic
%  Equation (DAE). 

%% Minimize time// energy consumption


% This code takes as inputs:

% Aerodynamic wind tunnel model 
% Propulsive  wind tunnel model
% Weights parameters from last version (June 2022)

%% Set references
global V1
V1 = 28 ; %Starting speed [m/s]
global V2 
V2 = 23 ; %Objective speed [m/s]
global D
D     = 0.8128; % Propeller Diameter [m]
global N_eng;
N_eng = 2; % Number of engines
global rho
rho   = 1.2133; % Air density [kg/m^3]

%% Aerodynamic model
p_CL_ac = [0.582;3.345731;-0.642635;-2.187085;0.713766;0.377985;-0.2946314;-0.004890;0.043976];


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


% Full model w.r.t angle of attack
CL = @(alpha_var) p_CL_ac(9)*alpha_var.^8+p_CL_ac(8)*alpha_var.^7 + p_CL_ac(7)*alpha_var.^6 + p_CL_ac(6)*alpha_var.^5 + p_CL_ac(5)*alpha_var.^4 + p_CL_ac(4)*alpha_var.^3 ...
    + p_CL_ac(3)*alpha_var.^2 + p_CL_ac(2)*alpha_var + p_CL_ac(1);

CD = @(alpha_var) p_CD_ac(10)*alpha_var.^9 + p_CD_ac(9)*alpha_var.^8 + p_CD_ac(8)*alpha_var.^7 + p_CD_ac(7)*alpha_var.^6 + p_CD_ac(6)*alpha_var.^5 + p_CD_ac(5)*alpha_var.^4 + p_CD_ac(4)*alpha_var.^3 ...
    + p_CD_ac(3)*alpha_var.^2 + p_CD_ac(2)*alpha_var + p_CD_ac(1);



% Full model w.r.t angle of attack

alpha_varvec = linspace(0,pi/2,180);


CLv = CL(alpha_varvec);

CL_max_w1_CR  = max(CLv);
Safety_F      = 1.2;
CL_max_w1_CR_ope = CL_max_w1_CR/Safety_F^2;




%% Loads and weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.5008; % Reference surface [m^2]

RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = @(V,alpha) 0.5*rho*V^2*S_ref*CL(alpha);
alpha_0 = fzero(@(alpha) L(V1,alpha)- W,5*pi/180)
%% For plot purposes

flag_plot = 0;

%% MINIMISATION

n_opt = fminbnd(@(n) cruise_DAEs_E(n,alpha_0,flag_plot), 55, 65)

%% For plot purposes

flag_plot = 1;

cruise_DAEs_E(n_opt,alpha_0,flag_plot);

%% Solve system for rps

function J_obj = cruise_DAEs_E(n,alpha_0,flag_plot)
global D
global N_eng
global rho
global V1
%% Define interpolations and initial conditions
Tlin = 200; % Steps for time
tspan = linspace(0,200,Tlin);      % Interpolate time of integration

y0 = [V1 0 alpha_0]';               % Initial conditions. X(0) = 0, V(0) = V1,  alpha_0(0) = alpha0
TOL = 1e-10;    
% [y0, yp0] = decic(@fun,t0,y0,1,yp0,0)

M = [1 0 0 ; 0 1 0; 0 0 0]; % Mass Matrix (to solver DAE)
%% Definition of options
optionsolver = odeset('MassSingular','yes','Mass',M,'AbsTol',TOL,'RelTol',TOL,'Events',@eventsH);

[t,y] = ode15s(@(t,y)fun(t,y,n),tspan,y0,optionsolver);


%% For energy calculations
mbatt  = 3; %[kg]
e0     = 720e3; %[J/kg]
E_batt = e0*mbatt;  %[J] Total energy of the battery packs
xf     = 644378/2.5; % Full Cruise range [m]
tau    = 0.2;

% Power coefficient (for Power calc)
CP0 = 0.0261518307541734;
CP1 = 0.0473735972985378;
CP2 = -0.16267474946046;
CP3 = 0.0247028469343899;

%% Results!!
t_vec     = interp1(linspace(0,1,numel(t)),t,linspace(0,1,Tlin));           % For time



V_vec     = interp1(linspace(0,1,numel(y(:,1))),y(:,1),linspace(0,1,Tlin)); % For V 
E_vec     = E_batt - interp1(linspace(0,1,numel(y(:,2))),y(:,2),linspace(0,1,Tlin)); % For energy 
P_vec     =  (CP3*V_vec.^3./(n.^3*D^3)+CP2*V_vec.^2./(n.^2*D^2)+CP1*V_vec./(n.*D)+CP0).*rho*N_eng*n.^3*D^5;

x_vec(1)  = 0;

for ii = 2:length(V_vec)
x_vec(ii) = V_vec (ii-1)*(t_vec(ii)-t_vec(ii-1)) + x_vec(ii-1);
end


if flag_plot ==1
    

    
    figure(1)
    plot(t_vec,x_vec,'k','LineWidth',1.2)
    grid on
    xlabel('Cruise time [s] ')
    ylabel('Horizontal distance [m]')
    Title1 = strcat('Horizontal distance vs time');
    title(Title1)
    
    
    figure(2)
    plot(t_vec,V_vec,'k','LineWidth',1.2)
    grid on
    hold on
    xlabel('Cruise time [s] ')
    ylabel('Flight speed [m/s]')
    Title1 = strcat('Flight speed vs time');
    title(Title1)
    legend('Total Speed V','Location','best')
    
    figure(3)
    plot(t_vec,E_vec,'k','LineWidth',1.2)
    grid on
    xlabel('Cruise time [s] ')
    ylabel('Energy consumption [m]')
    Title1 = strcat('Energy consumption vs time');
    title(Title1)
    
    
    figure(4)
    plot(t_vec,P_vec,'k','LineWidth',1.2)
    grid on
    hold on
    %yline(Pmax,'--b','LineWidth',1.2)
    xlabel('Cruise time [s] ')
    ylabel('Power consumption [W]')
    Title1 = strcat('Power consumption vs time');
    title(Title1)

    
end

E_cruise = E_vec(end); % Total energy consumed in climb
x_cruise = x_vec(end); % Total distance travelled in climb
J_obj = -E_cruise % Value of the objective function
end



function out = fun(t,y,n)
% Some parameters are not required within this function, but they are
% included anyway. Will delete in a further version of the code.

global D
global N_eng
global rho



%% Aerodynamic model
p_CL_ac = [0.582;3.345731;-0.642635;-2.187085;0.713766;0.377985;-0.2946314;-0.004890;0.043976];


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


% App. model w.r.t angle of attack
CL = @(alpha_var) p_CL_ac(9)*alpha_var.^8+p_CL_ac(8)*alpha_var.^7 + p_CL_ac(7)*alpha_var.^6 + p_CL_ac(6)*alpha_var.^5 + p_CL_ac(5)*alpha_var.^4 + p_CL_ac(4)*alpha_var.^3 ...
    + p_CL_ac(3)*alpha_var.^2 + p_CL_ac(2)*alpha_var + p_CL_ac(1);

CD = @(alpha_var) p_CD_ac(10)*alpha_var.^9 + p_CD_ac(9)*alpha_var.^8 + p_CD_ac(8)*alpha_var.^7 + p_CD_ac(7)*alpha_var.^6 + p_CD_ac(6)*alpha_var.^5 + p_CD_ac(5)*alpha_var.^4 + p_CD_ac(4)*alpha_var.^3 ...
    + p_CD_ac(3)*alpha_var.^2 + p_CD_ac(2)*alpha_var + p_CD_ac(1);



% App. model w.r.t angle of attack


alpha_varvec = linspace(0,pi/2,180);


CLv = CL(alpha_varvec);

CL_max_w1_CR  = max(CLv);
Safety_F      = 1.2;
CL_max_w1_CR_ope = CL_max_w1_CR/Safety_F^2;
%% Propulsive Model

CT0 =  0.0735880531010883;
CT1 = -0.0311758018412727;
CT2 = -0.249744726429543; 

%% Propulsive Model

CP0 = 0.0261518307541734;
CP1 = 0.0473735972985378;
CP2 = -0.16267474946046;
CP3 = 0.0247028469343899;

%% Electric
tau = 0.2;
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency




%% Loads and weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.5008; % Reference surface [m^2]
V_min_ope  = sqrt(2*W/(rho*S_ref*CL_max_w1_CR_ope));  % [m/s] Min operative speed

%% alpha_var relationship (from T = D)

%% Thrust
CT = @(y,n) CT2*y(1)^2/(n^2*D^2)+CT1*y(1)/(n*D)+CT0;

T       = @(y,n) CT(y,n)*N_eng*rho*n^2*D^4;

%% Lift

L = @(y) 0.5*rho*y(1)^2*S_ref*CL(y(3));

%% Drag

Dr = @(y) 0.5*rho*y(1)^2*S_ref.*CD(y(3));

%% Power
CP = @(y,n) CP3*y(1)^3/(n^3*D^3)+CP2*y(1)^2/(n^2*D^2)+CP1*y(1)/(n*D)+CP0;

P = @(y,n) CP(y,n)*N_eng*rho*n^3*D^5;

%% Define system of equations

% y2 = V
% y3 = E
% y4 = alpha_var

out = zeros(3,1);
out(1) = (T(y,n) - Dr(y))/mTOW ; % Long forces equation
out(2) = P(y,n)/(eta_m);
out(3) = L(y)-W;

end


function [position,isterminal,direction] = eventsH(t,y)
global V2
position = y(1)-V2; % The value that we want to be zero
isterminal = 1;           % Halt integration
direction = 0;            % The zero can be approached from either direction
end













