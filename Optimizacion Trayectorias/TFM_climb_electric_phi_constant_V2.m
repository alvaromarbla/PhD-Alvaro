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


n_vec_analysis = [53 55 57 59 61 63 65 67 69 71 73 75 rps_max];

for kk = 1:length(n_vec_analysis)
n_val = n_vec_analysis(kk);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% MINIMISATION
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

J = @(X) X(1)./(n_val*D);
phi = @(X) X(3)+X(4);

CP = @(X)  a1_P + b1_P*J(X) + c1_P*J(X).^2 + d1_P*J(X).^3 + e1_P*J(X).^4 + ...
    phi(X).*(b2_P*J(X) + c2_P*J(X).^2 + d2_P*J(X).^3) +  phi(X).^2.*(b3_P*J(X) + c3_P*J(X).^2)...
    +  phi(X).^3.*b4_P*J(X);

P_f = @(X) CP(X)*N_eng*rho*n_val^3*D^5;


%% For energy calculations
mbatt  = 3; %[kg]
e0     = 720e3; %[J/kg]
E_batt = e0*mbatt;  %[J] Total energy of the battery packs
xf     = 644378/4; % Full Cruise range [m]
tau    = 0.2;
C_cost = E_batt*(1-tau)/(xf); % Cost for objective function
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency

%% For plot purposes

E_v_fmincon = @(X) P_f(X)*deltaH/(eta_m*X(1)*sin(X(2)))-C_cost*deltaH*cot(X(2));
nonlcon = @constrains_CLIMB;

if kk == 1
x0_f = [20 5*pi/180 5*pi/180 5*pi/180];
else
x0_f = Xsol_fmincon;
end
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];

options = optimoptions("fmincon","Algorithm","interior-point",  "EnableFeasibilityMode",true, "SubproblemAlgorithm","cg");
[Xsol_fmincon,f] = fmincon(E_v_fmincon,x0_f,A,b,Aeq,beq,lb,ub,nonlcon,options)

%% Define interpolations and initial conditions
Tlin = 200; % Steps for time
t_end = (deltaH-delta0)/(Xsol_fmincon(1)*sin(Xsol_fmincon(2)));
t_vec = linspace(0,t_end,Tlin);      % Interpolate time of integration



%X(1) = V
%X(2) = gamma
%X(3) = alpha
%X(4) = epsilon


V_vec = Xsol_fmincon(1)*ones(1,length(t_vec));
gamma_var = Xsol_fmincon(2);
alpha_var = Xsol_fmincon(3);
eps_eng_val = Xsol_fmincon(4);
x_vec = V_vec.*t_vec*cos(gamma_var);
h_vec = V_vec.*t_vec*sin(gamma_var)+delta0;
E_vec = E_batt - P_f(Xsol_fmincon)*t_vec/eta_m;
n_vec = n_val*ones(1,length(t_vec));

P_vec = P_f(Xsol_fmincon)*ones(1,length(t_vec));

flag_plot = 0;

gamma_plot(kk) = gamma_var;
V_plot(kk)    = Xsol_fmincon(1);
alpha_plot(kk) = alpha_var;
epsilon_plot(kk) = Xsol_fmincon(4);
E_plot(kk) = E_vec(end)/E_batt*100;
J_sol (kk) = E_v_fmincon(Xsol_fmincon);
Pmax_Eng = 6.7e3; % Max Power per engine [W]
Pmax = N_eng*Pmax_Eng;


if flag_plot == 1

    figure(1)
    plot(t_vec,gamma_var*180/pi*ones(length(t_vec),1),'k','LineWidth',1.2)
    grid on
    xlabel('Climb time [s] ')
    ylabel('Trajectory angle [º]')
    Title1 = strcat(['Trajectory Angle \gamma vs time. Trajectory angle = ', num2str(gamma_var*180/pi), 'deg  \phi =', num2str((eps_eng_val+alpha_var)*180/pi)],'deg');
    title(Title1)

    figure(2)
    plot(t_vec,x_vec,'k','LineWidth',1.2)
    grid on
    xlabel('Climb time [s] ')
    ylabel('Horizontal distance [m]')
    Title1 = strcat(['Horizontal distance vs time. Trajectory angle = ', num2str(gamma_var*180/pi),'deg \phi =', num2str((eps_eng_val+alpha_var)*180/pi)],'deg');
    title(Title1)

    figure(3)
    plot(t_vec,h_vec,'k','LineWidth',1.2)
    grid on
    xlabel('Climb time [s] ')
    ylabel('Vertical distance [m]')
    Title1 = strcat(['Vertical distance vs time. Trajectory angle = ', num2str(gamma_var*180/pi), ' deg, \phi =', num2str((eps_eng_val+alpha_var)*180/pi)],'deg');
    title(Title1)

    figure(4)
    plot(t_vec,V_vec,'k','LineWidth',1.2)
    grid on
    hold on
    plot(t_vec,V_vec.*cos(gamma_var),'--b','LineWidth',1.2)
    plot(t_vec,V_vec.*sin(gamma_var),'-.r','LineWidth',1.2)
    xlabel('Climb time [s] ')
    ylabel('Flight speed [m/s]')
    Title1 = strcat(['Flight speed vs time. Trajectory angle = ', num2str(gamma_var*180/pi), 'deg, \phi =', num2str((eps_eng_val+alpha_var)*180/pi)],'deg');
    title(Title1)
    legend('Total Speed V','Horizontal Speed Vh','Vertical Speed Vv','Location','best')

    figure(5)
    plot(t_vec,E_vec,'k','LineWidth',1.2)
    grid on
    xlabel('Climb time [s] ')
    ylabel('Energy consumption [m]')
    Title1 = strcat(['Energy consumption vs time. Trajectory angle = ', num2str(gamma_var*180/pi), 'deg, \phi =', num2str((eps_eng_val+alpha_var)*180/pi)],'deg');
    title(Title1)


    figure(6)
    plot(t_vec,P_vec,'k','LineWidth',1.2)
    grid on
    hold on
    %yline(Pmax,'--b','LineWidth',1.2)
    xlabel('Climb time [s] ')
    ylabel('Power consumption [W]')
    Title1 = strcat(['Power consumption vs time. Engine revolutions = ', num2str(gamma_var*180/pi), 'deg, \phi =', num2str((eps_eng_val+alpha_var)*180/pi)],'deg');
    title(Title1)

    figure(7)
    plot(x_vec,h_vec,'k','LineWidth',1.2)
    grid on
    xlabel('Horizontal range [m] ')
    ylabel('Vertical distance [m]')
    Title1 = strcat(['Full Trajectory. Engine revolutions = ', num2str(gamma_var*180/pi), 'deg, \phi =', num2str((eps_eng_val+alpha_var)*180/pi)],'deg');
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

    L_vec =  0.5*rho.*V_vec.^2*S_ref*CL(alpha_var);

    %% Drag

    Dr_vec =  0.5*rho.*V_vec.^2*S_ref*CD(alpha_var);

    % Thrust

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

    T = @(V,n,alpha,eps_eng) CT(V,n,alpha,eps_eng).*N_eng*rho*n.^2*D^4;

    %t_vec = t_vec;
    %x_vec = x_vec;
    %h_vec = h_vec;
    %V_vec = V_vec;
    %n_vec = n_vec;
    gamma_vec = gamma_var*ones(1,length(t_vec));
    alpha_vec = alpha_var*ones(1,length(t_vec));
    epsilon_vec = eps_eng_val*ones(1,length(t_vec));
    phi_vec = alpha_vec+epsilon_vec;
    %L_vec = L_vec;
    T_vec = T(Xsol_fmincon(1),n_val,alpha_var,eps_eng_val)*ones(1,length(t_vec)) ;
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

end



figure(100)
plot(n_vec_analysis,V_plot,'k','LineWidth',1.2)
grid on
xlabel('Engine running speed [rps] ')
ylabel('Vertical distance [m]')
Title1 = strcat('Sensitivity Analysis. n vs V');
title(Title1)

figure(101)
plot(n_vec_analysis,alpha_plot*180/pi,'b','LineWidth',1.2)
grid on
hold on
plot(n_vec_analysis,epsilon_plot*180/pi,'-.r','LineWidth',1.2)
plot(n_vec_analysis,(epsilon_plot+alpha_plot)*180/pi,'--o','LineWidth',1.2)
xlabel('Engine running speed [rps] ')
ylabel('Angle [deg]')
legend('Angle of attack','Engine tilt angle','Engine inclination w.r.t. airflow')
Title1 = strcat('Sensitivity Analysis. n vs Relevant angles');
title(Title1)


figure(102)
plot(n_vec_analysis,gamma_plot*180/pi,'k','LineWidth',1.2)
grid on
xlabel('Engine running speed [rps] ')
ylabel('Trajetory angle [deg]')
Title1 = strcat('Sensitivity Analysis. n vs \gamma');
title(Title1)

figure(103)
plot(n_vec_analysis,E_plot,'k','LineWidth',1.2)
grid on
xlabel('Engine running speed [rps] ')
ylabel('Remaining segment energy [%]')
Title1 = strcat('Sensitivity Analysis. n vs Energy');
title(Title1)

figure(104)
plot(n_vec_analysis,J_sol,'k','LineWidth',1.2)
grid on
xlabel('Engine running speed [rps] ')
ylabel('Objective Function J [J]')
Title1 = strcat('Sensitivity Analysis. n vs J_{obj}');
title(Title1)


function [c,ceq] = constrains_CLIMB(X)
%% Propulsive model
N_eng = 2; % Number of engines
D     = 0.8124; % Propeller Diameter [m]
RPMMAX_APC = 150000; % Max RPM of AXI motors
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

n_val = 55;

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;

%% Loads and weights

mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S = 0.5008; % Reference surface [m^2]
rho   = 1.2133;

% Thrust

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

T = @(V,n,alpha,eps_eng) CT(V,n,alpha,eps_eng)*N_eng*rho*n^2*D^4;


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

%% Drag model
%% Aerodynamic model
p_CL_ac = [0.582;3.345731;-0.642635;-2.187085;0.713766;0.377985;-0.2946314;-0.004890;0.043976];


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


% Full model w.r.t angle of attack
CL = @(alpha) p_CL_ac(9)*alpha.^8+p_CL_ac(8)*alpha.^7 + p_CL_ac(7)*alpha.^6 + p_CL_ac(6)*alpha.^5 + p_CL_ac(5)*alpha.^4 + p_CL_ac(4)*alpha.^3 ...
    + p_CL_ac(3)*alpha.^2 + p_CL_ac(2)*alpha + p_CL_ac(1);

CD = @(alpha) p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + p_CD_ac(4)*alpha.^3 ...
    + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);

% Full model w.r.t angle of attack

alpha_varvec = linspace(0,pi/2,180);


CLv = CL(alpha_varvec);
CL_max_w1_CR  = max(CLv);

Safety_F      = 1.2;
CL_max_w1_CR_ope = CL_max_w1_CR/Safety_F^2;

V_min_ope  = sqrt(2*W/(rho*S*CL_max_w1_CR_ope));  % [m/s] Min operative speed

%X(1) = V
%X(2) = gamma
%X(3) = alpha
%X(4) = epsilon

% [Stall speed condition; Max RPS conditions]



c  = [-X(1)+V_min_ope;-X(2);X(2)-pi/2;T(X(1),n_val,X(3),X(4))-Tmax;phi(X(3),X(4))-pi/2;X(3)-pi/2;X(4)-pi/2;-CP(X(1),n_val,X(3),X(4))];

% T = D condition
ceq = [T(X(1),n_val,X(3),X(4)).*cos(phi(X(3),X(4)))-0.5*rho*X(1).^2*S*CD(X(3))-W*sin(X(2)); T(X(1),n_val,X(3),X(4)).*sin(phi(X(3),X(4)))+0.5*rho*X(1).^2*S*CL(X(3))-W*cos(X(2))];
end







