clc
clear all
close all

%% Calculate Heuristic Transition among 2 phases.
% Change criteria will be based on speed threshold



%% The first phase will be solved using DAEs.


% This code takes as inputs:

% Aerodynamic wind tunnel model
% Propulsive  wind tunnel model
% Weights parameters from last version (June 2022)

%% Model the Aircraft
[phase2_cond] = Heuristic_RevTrans;

function phase2_cond = Heuristic_RevTrans
%% Aerodynamic model

rho = 1.223;
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
W    = mTOW*g;
S_ref = 0.5008; % Reference surface [m^2]

%% For energy calculations
mbatt  = 3; %[kg]
e0     = 720e3; %[J/kg]
E_batt = e0*mbatt;  %[J] Total energy of the battery packs

V_min_ope  = sqrt(2*W/(rho*S_ref*CL_max_w1_CR_ope));  % [m/s] Min operative speed

N_eng = 2;  % Number of engines [-]
D     = 0.8128; % Propeller Diameter [m]
RPMMAX_APC = 150000; % Max RPM of AXI motors
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;

%% Flight conditions

delta0    = 100; % Initial altitude
t_estimate = 10; % Estimate of time, for calculating the phi_slope
eps_slope = pi/(2*t_estimate);
eps_0     = 0;
alpha_val = 0.2618; %0.10420 ; % Estimate of alpha for V = 23 m/s;0.133831601086466
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For plot purposes


[phase2_cond] = phase1_heuristicTrans(delta0, alpha_val,eps_slope, eps_0);

t_vec = phase2_cond(:,1);
x_vec = phase2_cond(:,2);
V_vec = phase2_cond(:,3);
E_vec = phase2_cond(:,4);
n_vec = phase2_cond(:,5);

%% Tilt rotor variation

eps_FF = pi/2;%6*pi/9 - alpha_val;
tF     = t_estimate;%(eps_FF - eps_0 )/eps_slope;
[~, t_index] = min(abs(t_vec-tF));
t_break = t_vec (t_index); %% Calculate phi break ()

eps_eng_vec = eps_0 + eps_slope.*t_vec+ (eps_FF - eps_0 - eps_slope.*t_vec).*heaviside(t_vec-tF);
alpha_vec   = ones(length(t_vec),1)*alpha_val;

J_vec   = V_vec./(n_vec*D); % Advance ratio
phi_vec = alpha_vec + eps_eng_vec;
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

CT_vec =  a1_T + b1_T*J_vec + c1_T*J_vec.^2 + d1_T*J_vec.^3 + e1_T*J_vec.^4 + ...
    phi_vec.*(b2_T*J_vec+ c2_T*J_vec.^2 + d2_T*J_vec.^3) +  phi_vec.^2.*(b3_T*J_vec + c3_T*J_vec.^2)...
    +  phi_vec.^3.*b4_T.*J_vec;

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

CP_vec =  a1_P + b1_P*J_vec + c1_P*J_vec.^2 + d1_P*J_vec.^3 + e1_P*J_vec.^4 + ...
    phi_vec.*(b2_P*J_vec+ c2_P*J_vec.^2 + d2_P*J_vec.^3) +  phi_vec.^2.*(b3_P*J_vec + c3_P*J_vec.^2)...
    +  phi_vec.^3.*b4_P.*J_vec;
%% Electric
tau   = 0.2;
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency

%% Loads and weight
mTOW  = 16.6; % Maximum T-O Mass [kg]
g     = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.5008; % Reference surface [m^2]


%% Thrust

T_vec   =  CT_vec.*N_eng*rho.*n_vec.^2*D^4;

%% Lift

CL_vec = CL(alpha_vec);
L_vec = 0.5*rho*V_vec.^2*S_ref.*CL_vec;

%% Drag
CD_vec = CD(alpha_vec);
Dr_vec =  0.5*rho*V_vec.^2*S_ref.*CD_vec;

%% Power

P_vec =  CP_vec.*N_eng*rho.*n_vec.^3*D^5;

figure(1)
plot(t_vec,x_vec,'LineWidth',1.5)
hold on
grid on
plot(t_break,x_vec(t_index),'or','LineWidth',1.5)
Title1 = strcat('Time [s] vs Distance [m]. ');
title(Title1)
xlabel('t [s]')
ylabel('x [m]')
legend('Horizontal distance profile','Engine Tilt breakpoint','Location','northwest')

figure(2)
plot(t_vec,V_vec,'LineWidth',1.5)
hold on
grid on
plot(t_break,V_vec(t_index),'or','LineWidth',1.5)
Title1 = strcat('Time [s] vs Speed [m/s]');
title(Title1)
xlabel('t [s]')
ylabel('V [m/s]')
legend('Velocity profile','Engine Tilt breakpoint','Location','northwest')

hFIG3 = figure(3)
fname = "Rev_trans_Eperc";
plot(t_vec,E_vec/E_batt*100,'LineWidth',1.5)
hold on
grid on
plot(t_break,E_vec(t_index)/E_batt*100,'or','LineWidth',1.5)
Title1 = strcat('Time [s] vs Energy perc [\%].');
title(Title1)
xlabel('t [s]')
ylabel('E perc [\%]')
legend('Energy consumption profile','Engine Tilt breakpoint','Location','northwest')

picturewidth = 20;
hw_ratio = 0.65;
set(findall(hFIG3,'-property','FontSize'),'FontSize',16)
set(findall(hFIG3,'-property','Box'),'Box','off') % optional
set(findall(hFIG3,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hFIG3,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hFIG3,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hFIG3,'Position');
set(hFIG3,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hFIG3,fname,'-dpdf','-vector','-fillpage')

figure(4)
plot(t_vec,n_vec,'LineWidth',1.5)
hold on
grid on
plot(t_break,n_vec(t_index),'or','LineWidth',1.5)
yline(rps_max,'--k','LineWidth',1.5)
Title1 = strcat('Time [s] vs Eng. revolutions [rps].');
title(Title1)
xlabel('t [s]')
ylabel('n [rps]')
legend('Engine revolutions distribution','Engine Tilt breakpoint','Location','best')

% figure(5)
% plot(t_vec,alpha_vec*180/pi,'LineWidth',1.5)
% Title1 = strcat(['Time [s] vs Angle of attack [deg]. Phase break at ' num2str(alpha_vec(iindex)*180/pi) ' º.']);
% title(Title1)
% hold on
% plot (t_vec(iindex),alpha_vec(iindex)*180/pi,'or','LineWidth',1.5)
% xlabel('t [s]')
% ylabel('\alpha [deg]')
% legend('Angle of attack distribution','Break point','Location','northeast')

figure(5)
plot(x_vec,CL_vec,'LineWidth',1.5)
hold on
grid on
Title1 = strcat('Horizontal distance [m] vs Aerodynamic Coefficients [-].');
title(Title1)
plot(x_vec,CD_vec,'LineWidth',1.5)
plot(t_break,CL_vec(t_index),'or','LineWidth',1.5)
plot(t_break,CD_vec(t_index),'^r','LineWidth',1.5)
xlabel('x [m]')
ylabel('Aerodynamic coefficients [-]')
legend('Lift Coefficient distribution over distance','Drag Coefficient distribution over distance','Engine Tilt breakpoint for C_L','Engine Tilt breakpoint for C_D','Location','northwest')


figure(6)
plot(x_vec,V_vec,'LineWidth',1.5)
Title1 = strcat('Horizontal distance [m] vs Speed [m/s].' );
title(Title1)
hold on
grid on
plot(x_vec(t_index),V_vec(t_index),'or','LineWidth',1.5)
xlabel('x [m]')
ylabel('V [m/s]')
legend('Velocity profile','Engine Tilt breakpoint','Location','northwest')

figure(7)
plot(x_vec,L_vec,'LineWidth',1.5)
Title1 = strcat('Horizontal distance [m] vs Aerodynamic Forces [N].');
hold on
grid on
title(Title1)
plot(x_vec,Dr_vec,'LineWidth',1.5)
plot(x_vec(t_index),L_vec(t_index),'or','LineWidth',1.5)
plot(x_vec(t_index),Dr_vec(t_index),'^r','LineWidth',1.5)
xlabel('x [m]')
ylabel('Aerodynamic forces [N]')
legend('Lift distribution over distance','Drag distribution over distance','Engine Tilt breakpoint for Lift','Engine Tilt breakpoint for Drag','Location','northwest')

hFIG8 = figure(8)
fname = "Rev_trans_P";

plot(x_vec,P_vec/1000,'LineWidth',1.5)
hold on
grid on
plot(x_vec(t_index),P_vec(t_index)/1000,'or','LineWidth',1.5)
Title1 = strcat('Horizontal distance [m] vs Power [kW].');
title(Title1)
xlabel('x [m]')
ylabel('P [kW]')
legend('Power distribution over distance','Engine Tilt breakpoint','Location','north')

picturewidth = 20;
hw_ratio = 0.65;
set(findall(hFIG8,'-property','FontSize'),'FontSize',16)
set(findall(hFIG8,'-property','Box'),'Box','off') % optional
set(findall(hFIG8,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hFIG8,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hFIG8,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hFIG8,'Position');
set(hFIG8,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hFIG8,fname,'-dpdf','-vector','-fillpage')

figure(9)
plot(x_vec,T_vec,'LineWidth',1.5)
hold on
grid on
plot(x_vec(t_index),T_vec(t_index),'or','LineWidth',1.5)
yline(Tmax,'--k','LineWidth',1.5)
Title1 = strcat('Horizontal distance [m] vs Thrust [N]. ');
title(Title1)
xlabel('x [m]')
ylabel('T [N]')
legend('Thrust distribution over distance','Engine Tilt breakpoint','Location','northwest')

figure(10)
plot(x_vec,phi_vec*180/pi,'LineWidth',1.5)
hold on
grid on
Title1 = strcat('Horizontal distance [m] vs Engine tilt angle [deg].');
title(Title1)
plot(x_vec,alpha_vec*180/pi,'LineWidth',1.5)
plot(x_vec,eps_eng_vec*180/pi,'LineWidth',1.5)
plot(x_vec(t_index),phi_vec(t_index)*180/pi,'or','LineWidth',1.5)
plot(x_vec(t_index),alpha_vec(t_index)*180/pi,'^r','LineWidth',1.5)
plot(x_vec(t_index),eps_eng_vec(t_index)*180/pi,'*r','LineWidth',1.5)

xlabel('x [m]')
ylabel('Relevant performance angles [deg]')
legend('Engine inclination w.r.t. airflow','Angle of attack','Engine tilt angle','Engine Tilt breakpoint for \phi','Engine Tilt breakpoint for \alpha','Engine Tilt breakpoint for \epsilon','Location','northwest')

figure(11)
plot(t_vec,phi_vec*180/pi,'LineWidth',1.5)
hold on
grid on
Title1 = strcat('Time [s] vs Engine tilt angle [deg].');
title(Title1)
plot(t_vec,alpha_vec*180/pi,'LineWidth',1.5)
plot(t_vec,eps_eng_vec*180/pi,'LineWidth',1.5)
plot(t_vec(t_index),phi_vec(t_index)*180/pi,'or','LineWidth',1.5)
plot(t_vec(t_index),alpha_vec(t_index)*180/pi,'^r','LineWidth',1.5)
plot(t_vec(t_index),eps_eng_vec(t_index)*180/pi,'*r','LineWidth',1.5)

xlabel('t [s]')
ylabel('Relevant performance angles [deg]')
legend('Engine inclination w.r.t. airflow','Angle of attack','Engine tilt angle','Engine Tilt breakpoint for \phi','Engine Tilt breakpoint for \alpha','Engine Tilt breakpoint for \epsilon','Location','northwest')



%t_vec = t_vec;
%x_vec = x_vec;
h_vec = 50*ones(length(t_vec),1);
%V_vec = V_vec;
%n_vec = n_vec;
gamma_vec = zeros(length(t_vec),1);
%alpha_vec = alpha_vec;
epsilon_vec = eps_eng_vec;
%phi_vec = phi_vec;
%L_vec = L_vec;
%T_vec = T_vec;
D_vec = Dr_vec;
%P_vec = P_vec;
E_vec =  E_batt-E_vec;
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


Mission_9 = [t_vec x_vec  h_vec V_vec n_vec gamma_vec alpha_vec  epsilon_vec phi_vec L_vec T_vec D_vec P_vec E_vec E_vec_perc];
save Mission_9.mat Mission_9

    function [phase2_cond,eps_F] = phase1_heuristicTrans (delta0, alpha_val, eps_slope, eps_0)

        %% Define interpolations and initial conditions
        Tlin = 300; % Steps for time
        tspan = linspace(0,60,Tlin);      % Interpolate time of integration

        y0 = [0 19.28 0 55]';               % Initial conditions. X(0) = 0, V(0) = V1,  E(0) = 0, n(0) = nCruise
        TOL = 1e-10;

        M = [1 0 0 0 ; 0 1 0 0; 0 0 1 0; 0 0 0 0]; % Mass Matrix (to solver DAE)

        %% Definition of options
        optionsolver = odeset('MassSingular','yes','Mass',M,'AbsTol',TOL,'RelTol',TOL,'Events',@eventsV);

        [t,y] = ode15s(@(t,y)fun1(t,y,alpha_val,eps_slope, eps_0),tspan,y0,optionsolver);

        phase2_cond = [t,y];
    end





    function out = fun1(t,y,alpha_val,eps_slope, eps_0)

        %% Aerodynamic model
        p_CL_ac = [0.582;3.345731;-0.642635;-2.187085;0.713766;0.377985;-0.2946314;-0.004890;0.043976];


        p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


        % App. model w.r.t angle of attack
        CL = @(alpha_var) p_CL_ac(9)*alpha_var.^8+p_CL_ac(8)*alpha_var.^7 + p_CL_ac(7)*alpha_var.^6 + p_CL_ac(6)*alpha_var.^5 + p_CL_ac(5)*alpha_var.^4 + p_CL_ac(4)*alpha_var.^3 ...
            + p_CL_ac(3)*alpha_var.^2 + p_CL_ac(2)*alpha_var + p_CL_ac(1);

        CD = @(alpha_var) p_CD_ac(10)*alpha_var.^9 + p_CD_ac(9)*alpha_var.^8 + p_CD_ac(8)*alpha_var.^7 + p_CD_ac(7)*alpha_var.^6 + p_CD_ac(6)*alpha_var.^5 + p_CD_ac(5)*alpha_var.^4 + p_CD_ac(4)*alpha_var.^3 ...
            + p_CD_ac(3)*alpha_var.^2 + p_CD_ac(2)*alpha_var + p_CD_ac(1);

        %% Propulsive Model
        % Nonzero components of the CT matrix

        %% Tilt rotor variation
        eps_FF = 5.3*pi/9 - alpha_val;
        tF     = (eps_FF - eps_0 )/eps_slope;
        eps_eng = @(t) eps_0 + eps_slope*t + (eps_FF - eps_0 - eps_slope*t)*heaviside(t-tF);

        J   = @(y) y(2)./(y(4)*D); % Advance ratio
        phi = @(alpha_var, t) alpha_var + eps_eng(t);

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

        CT = @(y,alpha,t)  a1_T + b1_T*J(y) + c1_T*J(y).^2 + d1_T*J(y).^3 + e1_T*J(y).^4 + ...
            phi(alpha,t).*(b2_T*J(y) + c2_T*J(y).^2 + d2_T*J(y).^3) +  phi(alpha,t).^2.*(b3_T*J(y) + c3_T*J(y).^2)...
            +  phi(alpha,t).^3.*b4_T*J(y);


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

        CP = @(y,alpha,t)  a1_P + b1_P*J(y) + c1_P*J(y).^2 + d1_P*J(y).^3 + e1_P*J(y).^4 + ...
            phi(alpha,t).*(b2_P*J(y) + c2_P*J(y).^2 + d2_P*J(y).^3) +  phi(alpha,t).^2.*(b3_P*J(y) + c3_P*J(y).^2)...
            +  phi(alpha,t).^3.*b4_P*J(y);

        %% Electric
        tau   = 0.2;
        eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency

        %% Loads and weight
        mTOW  = 16.6; % Maximum T-O Mass [kg]
        g     = 9.81; % Gravity [m/s^2]
        W     = mTOW*g;

        S_ref = 0.5008; % Reference surface [m^2]



        %% Thrust

        T   = @(y,alpha,t) CT(y,alpha,t)*N_eng*rho*y(4).^2*D^4;

        %% Lift

        L = @(y) 0.5*rho*y(2)^2*S_ref*CL(alpha_val);

        %% Drag

        Dr = @(y) 0.5*rho*y(2)^2*S_ref*CD(alpha_val);

        %% Power

        P = @(y,alpha,t) CP(y,alpha,t)*N_eng*rho*y(4).^3*D^5;

        %% Define system of equations

        % y1 = x
        % y2 = V
        % y3 = E
        % y4 = n

        %L([0,19.64,0,71])/(mTOW) + T([0,19.64,0,71],0.2618,0)*sin(phi(0.2618,0))/(mTOW) - g

        out    = zeros(4,1);
        out(1) = y(2); %dxdt = V
        out(2) = (T(y,alpha_val,t)*cos(phi(alpha_val,t))/mTOW - Dr(y))/mTOW; % Long forces equation
        out(3) = P(y,alpha_val,t)/(eta_m);
        out(4) = L(y)/(mTOW) + T(y,alpha_val,t)*sin(phi(alpha_val,t))/(mTOW) - g; % Transversal forces equation = 0


    end





    function [position,isterminal,direction] = eventsV(t,y)
        position = y(2)-0; % The value that we want to be zero
        isterminal = 1;           % Halt integration
        direction = 0;            % The zero can be approached from either direction
    end
end


