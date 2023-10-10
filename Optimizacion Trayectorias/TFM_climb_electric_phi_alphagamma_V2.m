clc
clear all
close all


%% Calculate solution to the problem:

%% Minimize time// energy consumption

%% We start from initial configuration (end of VTOL)
% and we finish at end of climb (cruise conditions)

% We will define one parameter law to be constant (for example, gamma)
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
deltaH = 350; %[m]
global delta0
delta0 = 50; %[m]
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


%% Set the engine inclination angle:

eps_eng_val = 10*pi/180;
%% Alpha value 

alpha_cruise = 6.601*pi/180*1.05; % Cruise angle of attack x1.05
%% For plot purposes

flag_plot = 0;

%% MINIMISATION
gamma_opt = 15*pi/180;

%gamma_opt = fminbnd(@(gamma_var) climb_calculator_E(gamma_var,eps_eng_val,alpha_cruise,flag_plot), 5*pi/180, 20*pi/180)

%% For plot purposes

flag_plot = 1;

climb_calculator_E(gamma_opt,eps_eng_val,alpha_cruise,flag_plot);

%% Solve system for rps

function J_obj = climb_calculator_E(gamma_var,eps_eng_val,alpha_cruise,flag_plot)
global D
global N_eng
global rho
global delta0
%% Define interpolations and initial conditions
Tlin = 200; % Steps for time
tspan = linspace(0,200,Tlin);      % Interpolate time of integration

y0 = [0 delta0 24.1209]';               % Initial conditions. X(0) = 0, Y(0) = h1; Vh(0) = Vstall,  Vv(0) = 0
TOL = 1e-10;                      % Define tolerance


%% Definition of options
optionsolver = odeset('AbsTol',TOL,'RelTol',TOL,'Events',@eventsH);

[t,y] = ode45(@(t,y)fun(t,y,gamma_var,eps_eng_val,alpha_cruise),tspan, y0,optionsolver);
% Results!
t_vec     = interp1(linspace(0,1,numel(t)),t,linspace(0,1,Tlin));           % For time



x_vec     = interp1(linspace(0,1,numel(y(:,1))),y(:,1),linspace(0,1,Tlin)); % For range 
h_vec     = interp1(linspace(0,1,numel(y(:,2))),y(:,2),linspace(0,1,Tlin)); % For height 
V_vec     = interp1(linspace(0,1,numel(y(:,3))),y(:,3),linspace(0,1,Tlin)); % For V 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loads and weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.5008; % Reference surface [m^2]
%% For energy calculations
mbatt  = 3; %[kg]
e0     = 720e3; %[J/kg]
E_batt = e0*mbatt;  %[J] Total energy of the battery packs
xf     = 644378/2.5; % Full Cruise range [m]
tau    = 0.2;
C_cost = E_batt/(xf*(1-tau)); % Cost for objective function
%% Propulsive Model
% Nonzero components of the CT matrix

J   = @(V,n) V./(n*D);
phi = @(alpha, eps_eng) alpha + eps_eng;

a1_T = 0.0735880531010883;
b1_T = -0.0311758018412727;
c1_T = -0.249744726429543;
d1_T = 0.143084420694372;
e1_T = 0.0261032283758581;

b2_T = -0.0982459868751664;
c2_T = 0.20127470719351;
d2_T = -0.173738749783189;

b3_T = 0.156239779501715;
c3_T = 0.0368592048084175;

b4_T = -0.0478709034281346;

CT = @(V,n,alpha,eps_eng)  a1_T + b1_T*J(V,n) + c1_T*J(V,n).^2 + d1_T*J(V,n).^3 + e1_T*J(V,n).^4 + ... 
    phi(alpha,eps_eng).*(b2_T*J(V,n) + c2_T*J(V,n).^2 + d2_T*J(V,n).^3) +  phi(alpha,eps_eng).^2.*(b3_T*J(V,n) + c3_T*J(V,n).^2)...
    +  phi(alpha,eps_eng).^3.*b4_T*J(V,n);
    

%% Propulsive Model

a1_P = 0.0261518307541734;
b1_P = 0.0473735972985378;
c1_P = -0.16267474946046;
d1_P = 0.0247028469343899;
e1_P = 0.0306053713439883;

b2_P = -0.0762350484603968;
c2_P = 0.148580471912353;
d2_P = -0.0726017200715775;

b3_P = 0.0897273366920878;
c3_P = 0.0122602815262456;

b4_P = -0.0486029866039398;

CP = @(V,n,alpha,eps_eng)  a1_P + b1_P*J(V,n) + c1_P*J(V,n).^2 + d1_P*J(V,n).^3 + e1_P*J(V,n).^4 + ... 
    phi(alpha,eps_eng).*(b2_P*J(V,n) + c2_P*J(V,n).^2 + d2_P*J(V,n).^3) +  phi(alpha,eps_eng).^2.*(b3_P*J(V,n) + c3_P*J(V,n).^2)...
    +  phi(alpha,eps_eng).^3.*b4_P*J(V,n);


%% Aerodynamic model
p_CL_ac = [0.582;3.345731;-0.642635;-2.187085;0.713766;0.377985;-0.2946314;-0.004890;0.043976];

% Full model w.r.t angle of attack
CL = @(alpha) p_CL_ac(9)*alpha.^8+p_CL_ac(8)*alpha.^7 + p_CL_ac(7)*alpha.^6 + p_CL_ac(6)*alpha.^5 + p_CL_ac(5)*alpha.^4 + p_CL_ac(4)*alpha.^3 ...
    + p_CL_ac(3)*alpha.^2 + p_CL_ac(2)*alpha + p_CL_ac(1);

%% Lift

L = @(V) 0.5*rho*V.^2*S_ref*CL(alpha_cruise);


%% Thrust
T   = @(V,n,alpha,eps_eng) CT(V,n,alpha,eps_eng).*N_eng*rho.*n.^2*D^4;

n0 = 60;
n_iter0 = n0*ones(Tlin,1);
options = optimset('TolX',1e-5,'MaxIter',1000);

n_vec = fsolve(@(n) T(V_vec,n,alpha_cruise,eps_eng_val) - (mTOW*g*cos(gamma_var)-L(V_vec))./sin(phi(alpha_cruise,eps_eng_val)),n_iter0,options)';
% size(n_vec)
% size(V_vec)
% %E_vec     = E_batt - interp1(linspace(0,1,numel(y(:,5))),y(:,5),linspace(0,1,Tlin)); % For energy 
% size(CP(V_vec,n_vec,alpha_cruise,eps_eng_val))
% size(CT(V_vec,n_vec,alpha_cruise,eps_eng_val))
% size(J(V_vec,n_vec))
T_vec =  CT(V_vec,n_vec,alpha_cruise,eps_eng_val).*rho*N_eng.*n_vec.^2*D^4;

P_vec =  CP(V_vec,n_vec,alpha_cruise,eps_eng_val).*rho*N_eng.*n_vec.^3*D^5;

E_vec = zeros(1,Tlin);
for kk = 1:Tlin
E_vec(kk) = E_batt - trapz(P_vec(1:kk))*t_vec(kk)/Tlin;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% % Power coefficient (for Power calc)
% 
% a1_P = 0.0261518307541734;
% b1_P = 0.0473735972985378;
% c1_P = -0.16267474946046;
% d1_P = 0.0247028469343899;
% e1_P = 0.0306053713439883;
% 
% b2_P = -0.0762350484603968;
% c2_P = 0.148580471912353;
% d2_P = -0.0726017200715775;
% 
% b3_P = 0.0897273366920878;
% c3_P = 0.0122602815262456;
% 
% b4_P = -0.0486029866039398;
% 
% 
% J = @(V,n) V./(n*D);
% phi = @(alpha, eps_eng) alpha + eps_eng;
% 
% CP = @(y,n,alpha, eps_eng)  a1_P + b1_P*J(y,n) + c1_P*J(y,n).^2 + d1_P*J(y,n).^3 + e1_P*J(y,n).^4 + ... 
%     phi(alpha, eps_eng).*(b2_P*J(y,n) + c2_P*J(y,n).^2 + d2_P*J(y,n).^3) +  phi(alpha, eps_eng).^2.*(b3_P*J(y,n) + c3_P*J(y,n).^2)...
%     +  phi(alpha, eps_eng).^3.*b4_P*J(y,n);

Pmax_Eng = 6.7e3; % Max Power per engine [W]
Pmax = N_eng*Pmax_Eng;

%% Results!!


%n_vec     = interp1(linspace(0,1,numel(y(:,4))),y(:,4),linspace(0,1,Tlin)); % For n 



if flag_plot ==1
    
    figure(1)
    plot(t_vec,gamma_var*180/pi*ones(length(t_vec),1),'k','LineWidth',1.2)
    grid on
    xlabel('Climb time [s] ')
    ylabel('Trajectory angle [º]')
    Title1 = strcat(['Trajectory Angle \gamma vs time. Trajectory angle = ', num2str(gamma_var*180/pi), 'deg  \phi =', num2str((eps_eng_val+alpha_cruise)*180/pi)],'deg');
    title(Title1)
    
    figure(2)
    plot(t_vec,x_vec,'k','LineWidth',1.2)
    grid on
    xlabel('Climb time [s] ')
    ylabel('Horizontal distance [m]')
    Title1 = strcat(['Horizontal distance vs time. Trajectory angle = ', num2str(gamma_var*180/pi),'deg \phi =', num2str((eps_eng_val+alpha_cruise)*180/pi)],'deg');
    title(Title1)
    
    figure(3)
    plot(t_vec,h_vec,'k','LineWidth',1.2)
    grid on
    xlabel('Climb time [s] ')
    ylabel('Vertical distance [m]')
    Title1 = strcat(['Vertical distance vs time. Trajectory angle = ', num2str(gamma_var*180/pi), ' deg, \phi =', num2str((eps_eng_val+alpha_cruise)*180/pi)],'deg');
    title(Title1)
    
    figure(4)
    plot(t_vec,V_vec,'k','LineWidth',1.2)
    grid on
    hold on
    plot(t_vec,V_vec.*cos(gamma_var),'--b','LineWidth',1.2)
    plot(t_vec,V_vec.*sin(gamma_var),'-.r','LineWidth',1.2)
    xlabel('Climb time [s] ')
    ylabel('Flight speed [m/s]')
    Title1 = strcat(['Flight speed vs time. Trajectory angle = ', num2str(gamma_var*180/pi), 'deg, \phi =', num2str((eps_eng_val+alpha_cruise)*180/pi)],'deg');
    title(Title1)
    legend('Total Speed V','Horizontal Speed Vh','Vertical Speed Vv','Location','best')
    
    figure(5)
    plot(t_vec,E_vec,'k','LineWidth',1.2)
    grid on
    xlabel('Climb time [s] ')
    ylabel('Energy consumption [m]')
    Title1 = strcat(['Energy consumption vs time. Trajectory angle = ', num2str(gamma_var*180/pi), 'deg, \phi =', num2str((eps_eng_val+alpha_cruise)*180/pi)],'deg');
    title(Title1)
    
    
    figure(6)
    plot(t_vec,P_vec,'k','LineWidth',1.2)
    grid on
    hold on
    %yline(Pmax,'--b','LineWidth',1.2)
    xlabel('Climb time [s] ')
    ylabel('Power consumption [W]')
    Title1 = strcat(['Power consumption vs time. Engine revolutions = ', num2str(gamma_var*180/pi), 'deg, \phi =', num2str((eps_eng_val+alpha_cruise)*180/pi)],'deg');
    title(Title1)
    
    figure(7)
    plot(x_vec,h_vec,'k','LineWidth',1.2)
    grid on
    xlabel('Horizontal range [m] ')
    ylabel('Vertical distance [m]')
    Title1 = strcat(['Full Trajectory. Engine revolutions = ', num2str(gamma_var*180/pi), 'deg, \phi =', num2str((eps_eng_val+alpha_cruise)*180/pi)],'deg');
    title(Title1)
    
%% Aerodynamic model
p_CL_ac = [0.582;3.345731;-0.642635;-2.187085;0.713766;0.377985;-0.2946314;-0.004890;0.043976];


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


% Full model w.r.t angle of attack
CL = @(alpha) p_CL_ac(9)*alpha.^8+p_CL_ac(8)*alpha.^7 + p_CL_ac(7)*alpha.^6 + p_CL_ac(6)*alpha.^5 + p_CL_ac(5)*alpha.^4 + p_CL_ac(4)*alpha.^3 ...
    + p_CL_ac(3)*alpha.^2 + p_CL_ac(2)*alpha + p_CL_ac(1);

CD = @(alpha) p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + p_CD_ac(4)*alpha.^3 ...
    + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);


%% Lift

L_vec =  0.5*rho.*V_vec.^2*S_ref*CL(alpha_cruise);

%% Drag

Dr_vec =  0.5*rho.*V_vec.^2*S_ref*CD(alpha_cruise);

%t_vec = t_vec;
%x_vec = x_vec;
%h_vec = h_vec;
%V_vec = V_vec;
%n_vec = n_vec;
gamma_vec = gamma_var*ones(1,length(t_vec));
alpha_vec = alpha_cruise*ones(1,length(t_vec));
epsilon_vec = eps_eng_val*ones(1,length(t_vec));
phi_vec = alpha_vec+epsilon_vec;
%L_vec = L_vec;
%T_vec = T_vec;
D_vec = Dr_vec;
%P_vec = P_vec;
%E_vec =  E_vec;
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


Mission_6 = [t_vec; x_vec;  h_vec; V_vec; n_vec; gamma_vec; alpha_vec;  epsilon_vec; phi_vec; L_vec; T_vec; D_vec; P_vec; E_vec; E_vec_perc]';
save Mission_6.mat Mission_6
% % 
% Mission_8 = [t_vec; x_vec;  h_vec; V_vec; n_vec; gamma_vec; alpha_vec;  epsilon_vec; phi_vec; L_vec; T_vec; D_vec; P_vec; E_vec; E_vec_perc]';
% save Mission_8.mat Mission_8




end

E_climb = E_vec(end); % Total energy consumed in climb
x_climb  = x_vec(end); % Total distance travelled in climb
J_obj = -E_climb + x_climb*C_cost % Value of the objective function
end



function dy = fun(t,y,gamma_var,eps_eng_val,alpha_cruise)
% Some parameters are not required within this function, but they are
% included anyway. Will delete in a further version of the code.

global rho



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


%% Electric
tau = 0.2;
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency




%% Loads and weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.5008; % Reference surface [m^2]
V_min_ope  = sqrt(2*W/(rho*S_ref*CL_max_w1_CR_ope));  % [m/s] Min operative speed

%% Lift

L = @(y) 0.5*rho*y(3)^2*S_ref*CL(alpha_cruise);

%% Drag

Dr = @(y) 0.5*rho*y(3)^2*S_ref*CD(alpha_cruise);

%% Engine inclination angle
phi = @(alpha, eps_eng) alpha + eps_eng;

%% Define system of equations

% y1 = x
% y2 = h
% y3 = V


dy = zeros(3,1);
dy(1) = y(3)*cos(gamma_var); %dxdt = Vcosgamma
dy(2) = y(3)*sin(gamma_var); % dh/dt = Vsingamma
dy(3) = (g*cos(gamma_var)-L(y)/mTOW)*cot(phi(alpha_cruise,eps_eng_val)) - Dr(y)/mTOW -g*sin(gamma_var); % Long forces equation (V)
%dy(4) = L(y)/(mTOW*y(3)) + T(y,alpha_cruise,eps_eng_val)*sin(phi(alpha_cruise,eps_eng_val))/(mTOW*y(3)) - g*cos(gamma_var)/y(3); % Transversal forces equation (gamma);
%dy(5) = P(y,alpha_cruise,eps_eng_val)/((1-tau)*eta_m); %(E)

end


function [position,isterminal,direction] = eventsH(t,y)
global deltaH
position = y(2)-deltaH; % The value that we want to be zero
isterminal = 1;           % Halt integration
direction = +1;            % The zero can be approached from either direction
end







