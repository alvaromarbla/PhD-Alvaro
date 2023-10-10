clc
clear all
close all

%% Calculate cruise from V1 to V2 given a certain RPS. The angle of attack and speed will be time dependant.
%  This problem will be solved in an analytical manner using integration.
%  The solution requires the integration of a Differential Algebraic
%  Equation (DAE). This will be easier as the algebraic equation can be
%  expressed in an explicit manner in terms of one of the variables.

%% Minimize time// energy consumption


% This code takes as inputs:

% Aerodynamic wind tunnel model // polar simplified model from previous
% Propulsive  wind tunnel model
% Weights parameters from last version (June 2022)

%% Set references
global V1
V1 = 24.12 ; %Starting speed [m/s]
global V2 
V2 = 19.28; %Objective speed [m/s]
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


% Formulate aero drag polar

CL0      = p_CL_ac(1);
CL_alpha = p_CL_ac(2);

CD2      = p_CD_ac(3)/CL_alpha^2; % CLO
CD1      = p_CD_ac(2)/CL_alpha - 2*CD2*CL0; % K1
CD0      = p_CD_ac(1)-CD1*CL0-CD2*CL0^2; % K2


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

%n_opt = fminbnd(@(n) cruise_calculator_E(n,flag_plot), 60, 78)
n_opt = 52
%% For plot purposes

flag_plot = 1;

cruise_calculator_E(n_opt,flag_plot);

%% Solve system for rps

function J_obj = cruise_calculator_E(n,flag_plot)
global D
global N_eng
global rho
global V1
%% Define interpolations and initial conditions
Tlin = 200; % Steps for time
tspan = linspace(0,200,Tlin);      % Interpolate time of integration

y0 = [0 V1 0]';               % Initial conditions. X(0) = 0, V(0) = V1,  E(0) = 0
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

%% Loads and weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.5008; % Reference surface [m^2]
% Power coefficient (for Power calc)
CP0 = 0.0261518307541734;
CP1 = 0.0473735972985378;
CP2 = -0.16267474946046;
CP3 = 0.0247028469343899;

%% Results!!
t_vec     = interp1(linspace(0,1,numel(t)),t,linspace(0,1,Tlin));           % For time



x_vec     = interp1(linspace(0,1,numel(y(:,1))),y(:,1),linspace(0,1,Tlin)); % For range 
V_vec     = interp1(linspace(0,1,numel(y(:,2))),y(:,2),linspace(0,1,Tlin)); % For V 
E_vec     = E_batt - interp1(linspace(0,1,numel(y(:,3))),y(:,3),linspace(0,1,Tlin)); % For energy 
P_vec     =  (CP3*V_vec.^3./(n.^3*D^3)+CP2*V_vec.^2./(n.^2*D^2)+CP1*V_vec./(n.*D)+CP0).*rho*N_eng*n.^3*D^5;

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

%% Calc L, D
%% Aerodynamic model
p_CL_ac = [0.582;3.345731;-0.642635;-2.187085;0.713766;0.377985;-0.2946314;-0.004890;0.043976];


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];

% App. model w.r.t angle of attack
CL = @(alpha) p_CL_ac(2)*alpha + p_CL_ac(1);

CD = @(alpha) p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);

% Formulate aero drag polar

CL0      = p_CL_ac(1);
CL_alpha = p_CL_ac(2);

CD2      = p_CD_ac(3)/CL_alpha^2; % CLO
CD1      = p_CD_ac(2)/CL_alpha - 2*CD2*CL0; % K1
CD0      = p_CD_ac(1)-CD1*CL0-CD2*CL0^2; % K2


alpha_vec =  (2*W./(rho*V_vec.^2*S_ref)-CL0)/CL_alpha;

L_vec = 0.5*rho.*V_vec.^2.*S_ref.*CL(alpha_vec);


Dr_vec = 0.5*rho.*V_vec.^2.*S_ref.*CD(alpha_vec);



%% Propulsive Model

CT0 =  0.0735880531010883;
CT1 = -0.0311758018412727;
CT2 = -0.249744726429543; 

%% Thrust
CT = @(y,n) CT2*y.^2./(n.^2*D^2)+CT1.*y./(n*D)+CT0;

T_vec   = CT(V_vec,n).*N_eng*rho*n.^2*D^4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%t_vec = t_vec;
%x_vec = x_vec;
h_vec = 50*ones(1,length(t_vec));
%V_vec = V_vec;
n_vec = n*ones(1,length(t_vec));
gamma_vec = zeros(1,length(t_vec));
%alpha_vec = alpha_vec;
epsilon_vec = zeros(1,length(t_vec));
phi_vec = alpha_vec;
%L_vec = L_vec;
%T_vec = T_vec;
D_vec = Dr_vec;
%P_vec = P_vec;
%E_vec = E_vec;
E_vec_perc = (E_vec)/E_batt*100;
% 
% size(t_vec)
% size(x_vec)
% size(h_vec)
% size(V_vec)
% size(n_vec)
% size(gamma_vec)
% size(alpha_vec)
% size(epsilon_vec)
% size(phi_vec)
% size(L_vec)
% size(D_vec)
% size(T_vec)
% size(P_vec)
% size(E_vec)
% size(E_vec_perc)


Mission_5c = [t_vec; x_vec ; h_vec; V_vec; n_vec; gamma_vec; alpha_vec ; epsilon_vec; phi_vec; L_vec; T_vec; D_vec; P_vec; E_vec; E_vec_perc]';
save Mission_5c.mat Mission_5c


    
end

E_cruise = E_vec(end); % Total energy consumed in climb
x_cruise = x_vec(end); % Total distance travelled in climb
J_obj = -E_cruise % Value of the objective function
end



function dy = fun(t,y,n)
% Some parameters are not required within this function, but they are
% included anyway. Will delete in a further version of the code.

global D
global N_eng
global rho



%% Aerodynamic model
p_CL_ac = [0.582;3.345731;-0.642635;-2.187085;0.713766;0.377985;-0.2946314;-0.004890;0.043976];


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


% App. model w.r.t angle of attack
CL = @(alpha) p_CL_ac(2)*alpha + p_CL_ac(1);

CD = @(alpha) p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);


% Formulate aero drag polar

CL0      = p_CL_ac(1);
CL_alpha = p_CL_ac(2);

CD2      = p_CD_ac(3)/CL_alpha^2; % CLO
CD1      = p_CD_ac(2)/CL_alpha - 2*CD2*CL0; % K1
CD0      = p_CD_ac(1)-CD1*CL0-CD2*CL0^2; % K2

% App. model w.r.t angle of attack


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

%% Alpha relationship (from T = D)

alpha_fun = @(y) (2*W./(rho*y(2).^2*S_ref)-CL0)/CL_alpha;

%% Thrust
CT = @(y,n) CT2*y(2)^2/(n^2*D^2)+CT1*y(2)/(n*D)+CT0;

T       = @(y,n) CT(y,n)*N_eng*rho*n^2*D^4;

% %% Lift
% 
% L = @(y) 0.5*rho*y(3)^2*S_ref*CL(alpha_fun(y));

%% Drag

Dr = @(y) 0.5*rho*y(2)^2*S_ref.*CD(alpha_fun(y));

%% Power
CP = @(y,n) CP3*y(2)^3/(n^3*D^3)+CP2*y(2)^2/(n^2*D^2)+CP1*y(2)/(n*D)+CP0;

P = @(y,n) CP(y,n)*N_eng*rho*n^3*D^5;

%% Define system of equations

% y1 = x
% y2 = V
% y3 = E


dy = zeros(3,1);
dy(1) = y(2); %dxdt = V
dy(2) = (T(y,n) - Dr(y))/mTOW ; % Long forces equation
dy(3) = P(y,n)/(eta_m);
end


function [position,isterminal,direction] = eventsH(t,y)
global V2
position = y(2)-V2; % The value that we want to be zero
isterminal = 1;           % Halt integration
direction = -1;            % The zero can be approached from either direction
end













