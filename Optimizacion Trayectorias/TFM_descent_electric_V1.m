clc
clear all
close all

%% Calculate solution to the problem:

%% Minimize time// energy consumption

%% We start from initial configuration (end of VTOL)
% and we finish at end of climb (cruise conditions)

% We will define one parameter law to be constant (for example, n)
% and check the evolution of the other variables as a function of time.

% For this first approach alpha is set to constant, phi = 0.

% We will need to integrate the function dynamics, and the resulting energy
% /time will arise as a byproduct of this method.

% This code takes as inputs:

% Aerodynamic wind tunnel model // polar simplified model from previous
% Propulsive  wind tunnel model
% Weights parameters from last version (June 2022)

%% Set deltaH (altitude difference)
global deltaH
deltaH = 10; %[m]
global delta0
delta0 = 310; %[m]
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
CL = @(alpha) p_CL_ac(9)*alpha.^8+p_CL_ac(8)*alpha.^7 + p_CL_ac(7)*alpha.^6 + p_CL_ac(6)*alpha.^5 + p_CL_ac(5)*alpha.^4 + p_CL_ac(4)*alpha.^3 ...
    + p_CL_ac(3)*alpha.^2 + p_CL_ac(2)*alpha + p_CL_ac(1);

CD = @(alpha) p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + p_CD_ac(4)*alpha.^3 ...
    + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);

% Full model w.r.t angle of attack

alphavec = linspace(0,pi/2,180);


CLv = CL(alphavec);

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




%% For plot purposes

flag_plot = 0;

%% MINIMISATION

%n_opt = fminbnd(@(n) climb_calculator_E(n,flag_plot), 55, 75)
n_opt = 50;
%% For plot purposes

flag_plot = 1;

climb_calculator_E(n_opt,flag_plot);

%% Solve system for rps

function J_obj = climb_calculator_E(n,flag_plot)
global D
global N_eng
global rho
global delta0
%% Define interpolations and initial conditions
Tlin = 200; % Steps for time
tspan = linspace(0,200,Tlin);      % Interpolate time of integration

y0 = [0 delta0 21 -0.01 0]';               % Initial conditions. X(0) = 0, Y(0) = h1; Vh(0) = Vstall,  Vv(0) = 0
TOL = 1e-10;                      % Define tolerance

%% Definition of options
optionsolver = odeset('AbsTol',TOL,'RelTol',TOL,'Events',@eventsH);

[t,y] = ode45(@(t,y)fun(t,y,n),tspan,y0,optionsolver);

%% For energy calculations
mbatt  = 3; %[kg]
e0     = 720e3; %[J/kg]
E_batt = e0*mbatt;  %[J] Total energy of the battery packs
xf     = 644378/2.5; % Full Cruise range [m]
tau    = 0.2;
C_cost = E_batt/(xf*(1-tau)); % Cost for objective function

% Power coefficient (for Power calc)

CP0 = 0.0261518307541734;
CP1 = 0.0473735972985378;
CP2 = -0.16267474946046;
CP3 = 0.0247028469343899;

Pmax_Eng = 6.7e3; % Max Power per engine [W]
Pmax = N_eng*Pmax_Eng;

%% Results!!

t_vec     = interp1(linspace(0,1,numel(t)),t,linspace(0,1,Tlin));           % For time



x_vec     = interp1(linspace(0,1,numel(y(:,1))),y(:,1),linspace(0,1,Tlin)); % For range 
h_vec     = interp1(linspace(0,1,numel(y(:,2))),y(:,2),linspace(0,1,Tlin)); % For height 
V_vec     = interp1(linspace(0,1,numel(y(:,3))),y(:,3),linspace(0,1,Tlin)); % For V 
gamma_vec = interp1(linspace(0,1,numel(y(:,4))),y(:,4),linspace(0,1,Tlin)); % For gamma 
E_vec     = E_batt - interp1(linspace(0,1,numel(y(:,5))),y(:,5),linspace(0,1,Tlin)); % For energy 
P_vec     =  (CP3*V_vec.^3./(n.^3*D^3)+CP2*V_vec.^2./(n.^2*D^2)+CP1*V_vec./(n.*D)+CP0).*rho*N_eng*n.^3*D^5;


if flag_plot ==1
    
    figure(1)
    plot(t_vec,gamma_vec*180/pi,'k','LineWidth',1.2)
    xlabel('Climb time [s] ')
    ylabel('Trajectory angle [º]')
    Title1 = strcat(['Trajectory Angle \gamma vs time. Engine revolutions = ', num2str(n), 'rps']);
    title(Title1)
    
    figure(2)
    plot(t_vec,x_vec,'k','LineWidth',1.2)
    xlabel('Climb time [s] ')
    ylabel('Horizontal distance [m]')
    Title1 = strcat(['Horizontal distance vs time. Engine revolutions = ', num2str(n),'rps']);
    title(Title1)
    
    figure(3)
    plot(t_vec,h_vec,'k','LineWidth',1.2)
    xlabel('Climb time [s] ')
    ylabel('Vertical distance [m]')
    Title1 = strcat(['Vertical distance vs time. Engine revolutions = ', num2str(n), 'rps']);
    title(Title1)
    
    figure(4)
    plot(t_vec,V_vec,'k','LineWidth',1.2)
    hold on
    plot(t_vec,V_vec.*cos(gamma_vec),'--b','LineWidth',1.2)
    plot(t_vec,V_vec.*sin(gamma_vec),'-.r','LineWidth',1.2)
    xlabel('Climb time [s] ')
    ylabel('Flight speed [m/s]')
    Title1 = strcat(['Flight speed vs time. Engine revolutions = ', num2str(n), 'rps']);
    title(Title1)
    legend('Total Speed V','Horizontal Speed Vh','Vertical Speed Vv','Location','best')
    
    figure(5)
    plot(t_vec,E_vec,'k','LineWidth',1.2)
    xlabel('Climb time [s] ')
    ylabel('Energy consumption [m]')
    Title1 = strcat(['Energy consumption vs time. Engine revolutions = ', num2str(n), 'rps']);
    title(Title1)
    
    
    figure(6)
    plot(t_vec,P_vec,'k','LineWidth',1.2)
    hold on
    %yline(Pmax,'--b','LineWidth',1.2)
    xlabel('Climb time [s] ')
    ylabel('Power consumption [W]')
    Title1 = strcat(['Power consumption vs time. Engine revolutions = ', num2str(n), 'rps']);
    title(Title1)
    
    figure(7)
    plot(x_vec,h_vec,'k','LineWidth',1.2)
    xlabel('Horizontal range [m] ')
    ylabel('Vertical distance [m]')
    Title1 = strcat(['Full Trajectory. Engine revolutions = ', num2str(n), 'rps']);
    title(Title1)
    
end

E_climb = E_vec(end); % Total energy consumed in climb
x_climb  = x_vec(end); % Total distance travelled in climb
J_obj = -E_climb + x_climb*C_cost % Value of the objective function
end



function dy = fun(t,y,n)
% Some parameters are not required within this function, but they are
% included anyway. Will delete in a further version of the code.

global D
global N_eng
global rho
%% Aerodynamic model

alpha_cruise = 6.601*pi/180*1.0; % Cruise angle of attack x1.05


%% Aerodynamic model
p_CL_ac = [0.582;3.345731;-0.642635;-2.187085;0.713766;0.377985;-0.2946314;-0.004890;0.043976];


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


% Full model w.r.t angle of attack
CL = @(alpha) p_CL_ac(9)*alpha.^8+p_CL_ac(8)*alpha.^7 + p_CL_ac(7)*alpha.^6 + p_CL_ac(6)*alpha.^5 + p_CL_ac(5)*alpha.^4 + p_CL_ac(4)*alpha.^3 ...
    + p_CL_ac(3)*alpha.^2 + p_CL_ac(2)*alpha + p_CL_ac(1);

CD = @(alpha) p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + p_CD_ac(4)*alpha.^3 ...
    + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);
% Full model w.r.t angle of attack


alphavec = linspace(0,pi/2,180);


CLv = CL(alphavec);

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


%% Thrust
CT = @(y,n) CT2*y(3)^2/(n^2*D^2)+CT1*y(3)/(n*D)+CT0;

T       = @(y,n) CT(y,n)*N_eng*rho*n^2*D^4;

%% Lift

L = @(y) 0.5*rho*y(3)^2*S_ref*CL(alpha_cruise);

%% Drag

Dr = @(y) 0.5*rho*y(3)^2*S_ref*CD(alpha_cruise);

%% Power
CP = @(y,n) CP3*y(3)^3/(n^3*D^3)+CP2*y(3)^2/(n^2*D^2)+CP1*y(3)/(n*D)+CP0;

P = @(y,n) CP(y,n)*N_eng*rho*n^3*D^5;

%% Define system of equations

% y1 = x
% y2 = h
% y3 = V
% y4 = gamma
% y5 = E


dy = zeros(5,1);
dy(1) = y(3)*cos(y(4)); %dxdt = Vcosgamma
dy(2) = y(3)*sin(y(4)); % dh/dt = Vsingamma
dy(3) = (T(y,n) - Dr(y))/mTOW -g*sin(y(4)); % Long forces equation
dy(4) = L(y)/(mTOW*y(3)) - g*cos(y(4))/y(3); % Transversal forces equation
dy(5) = P(y,n)/((1-tau)*eta_m);
end


function [position,isterminal,direction] = eventsH(t,y)
global deltaH
position = y(2)-deltaH; % The value that we want to be zero
isterminal = 1;           % Halt integration
direction = 0;            % The zero can be approached from either direction
end







