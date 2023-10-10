%% FULL MISSION !!

% Gathers all data in .mat format from previous missions, creates full
% matrix and represents the data 


clear all
close all
clc

%% Legend for variables:

% Mission_1 = t
% Mission_2 = x
% Mission_3 = h
% Mission_4 = V
% Mission_5 = n
% Mission_6 = gamma
% Mission_7 = alpha
% Mission_8 = epsilon
% Mission_9 = phi
% Mission_10 = L
% Mission_11 = T
% Mission_12 = D
% Mission_13 = P
% Mission_14 = E
% Mission_15 = E_%





%% Legend for missions:

% 1 : Taxi (T = 0.6*W) 30 s
% 2 : VTO
% 3 : Hover 30 s
% 4 : V to H Trans
% 5 : Acelerate till cruise speed
% 6 : Plane mode climb gamma cte
% 5b : Acelerate till cruise speed
% 7 : Cruise extended range
% 5c : Acelerate till cruise speed
% 8 : Descent gamma cte
% 9 : H to V Trans
% 10 : Hover 30 s
% 11 : V Land

%% Electric
tau = 0.0; 
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency
mbatt = 3; %[kg]
e0    = 720e3; %[J/kg]
%E     = e0*mbatt;  %[J] Total energy of the battery packs
E = 1.170126925152786e+06 - e0*mbatt*0.2; % Total energy after mission to perform extended cruise
N_eng = 2; % Number of engines
D     = 0.8128; % Propeller Diameter [m]

RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;


% Load missions

%% Taxi
load('Mission_1.mat')
t_1 = Mission_1(:,1);
x_1 = Mission_1(:,2);
h_1 = Mission_1(:,3);
V_1 = Mission_1(:,4);
n_1 = Mission_1(:,5);
gamma_1 = Mission_1(:,6);
alpha_1 = Mission_1(:,7);
epsilon_1 = Mission_1(:,8);
phi_1 = Mission_1(:,9);
L_1 = Mission_1(:,10);
T_1 = Mission_1(:,11);
D_1 = Mission_1(:,12);
P_1 = Mission_1(:,13);
E_1 = Mission_1(:,14); 
Eperc_1 = Mission_1(:,15);

E_ini = e0*mbatt;
%% VTO 

load('Mission_2.mat')
t_2 = Mission_2(:,1)+t_1(end);
x_2 = Mission_2(:,2)+x_1(end);
h_2 = Mission_2(:,3);
V_2 = Mission_2(:,4);
n_2 = Mission_2(:,5);
gamma_2 = Mission_2(:,6);
alpha_2 = Mission_2(:,7);
epsilon_2 = Mission_2(:,8);
phi_2 = Mission_2(:,9);
L_2 = Mission_2(:,10);
T_2 = Mission_2(:,11);
D_2 = Mission_2(:,12);
P_2 = Mission_2(:,13);
E_2 = Mission_2(:,14) - (E_ini-E_1(end)) ; 
Eperc_2 = Mission_2(:,15) -(100 -Eperc_1(end)) ;


%% Hover

load('Mission_3.mat')
t_3 = Mission_3(:,1)+t_2(end);
x_3 = Mission_3(:,2)+x_2(end);
h_3 = Mission_3(:,3)+50;
V_3 = Mission_3(:,4);
n_3 = Mission_3(:,5);
gamma_3 = Mission_3(:,6);
alpha_3 = Mission_3(:,7);
epsilon_3 = Mission_3(:,8);
phi_3 = Mission_3(:,9);
L_3 = Mission_3(:,10);
T_3 = Mission_3(:,11);
D_3 = Mission_3(:,12);
P_3 = Mission_3(:,13);
E_3 = Mission_3(:,14) - (E_ini-E_2(end)) ; 
Eperc_3 = Mission_3(:,15) -(100 -Eperc_2(end)) ;


%% V to H trans
load('Mission_4.mat')
t_4 = Mission_4(:,1)+t_3(end);
x_4 = Mission_4(:,2)+x_3(end);
h_4 = Mission_4(:,3);
V_4 = Mission_4(:,4);
n_4 = Mission_4(:,5);
gamma_4 = Mission_4(:,6);
alpha_4 = Mission_4(:,7);
epsilon_4 = Mission_4(:,8);
phi_4 = Mission_4(:,9);
L_4 = Mission_4(:,10);
T_4 = Mission_4(:,11);
D_4 = Mission_4(:,12);
P_4 = Mission_4(:,13);
E_4 = Mission_4(:,14) - (E_ini-E_3(end)) ; 
Eperc_4 = Mission_4(:,15) -(100 -Eperc_3(end)) ;

%% Acceleration

load('Mission_5.mat')
t_5 = Mission_5(:,1)+t_4(end); %%
x_5 = Mission_5(:,2)+x_4(end);
h_5 = Mission_5(:,3);
V_5 = Mission_5(:,4);
n_5 = Mission_5(:,5);
gamma_5 = Mission_5(:,6);
alpha_5 = Mission_5(:,7);
epsilon_5 = Mission_5(:,8);
phi_5 = Mission_5(:,9);
L_5 = Mission_5(:,10);
T_5 = Mission_5(:,11);
D_5 = Mission_5(:,12);
P_5 = Mission_5(:,13);
E_5 = Mission_5(:,14) - (E_ini-E_4(end)) ;  %%
Eperc_5 = Mission_5(:,15) -(100 -Eperc_4(end)) ; %%

%% Climb

load('Mission_6.mat')
t_6 = Mission_6(:,1)+t_5(end); %%
x_6 = Mission_6(:,2)+x_5(end);
h_6 = Mission_6(:,3);
V_6 = Mission_6(:,4);
n_6 = Mission_6(:,5);
gamma_6 = Mission_6(:,6);
alpha_6 = Mission_6(:,7);
epsilon_6 = Mission_6(:,8);
phi_6 = Mission_6(:,9);
L_6 = Mission_6(:,10);
T_6 = Mission_6(:,11);
D_6 = Mission_6(:,12);
P_6 = Mission_6(:,13);
E_6 = Mission_6(:,14) - (E_ini-E_5(end)) ;  %%
Eperc_6 = Mission_6(:,15) -(100 -Eperc_5(end)) ; %%

%% Acceleration 2

load('Mission_5b.mat')
t_5b = Mission_5b(:,1)+t_6(end); %%
x_5b = Mission_5b(:,2)+x_6(end);
h_5b = Mission_5b(:,3)+300;
V_5b = Mission_5b(:,4);
n_5b = Mission_5b(:,5);
gamma_5b = Mission_5b(:,6);
alpha_5b = Mission_5b(:,7);
epsilon_5b = Mission_5b(:,8);
phi_5b = Mission_5b(:,9);
L_5b = Mission_5b(:,10);
T_5b = Mission_5b(:,11);
D_5b = Mission_5b(:,12);
P_5b = Mission_5b(:,13);
E_5b = Mission_5b(:,14) - (E_ini-E_6(end)) ;  %%
Eperc_5b = Mission_5b(:,15) -(100 -Eperc_6(end)) ; %%

%% Cruise
load('Mission_7.mat')
t_7 = Mission_7(:,1)+t_5b(end); %%
x_7 = Mission_7(:,2)+x_5b(end);
h_7 = Mission_7(:,3);
V_7 = Mission_7(:,4);
n_7 = Mission_7(:,5);
gamma_7 = Mission_7(:,6);
alpha_7 = Mission_7(:,7);
epsilon_7 = Mission_7(:,8);
phi_7 = Mission_7(:,9);
L_7 = Mission_7(:,10);
T_7 = Mission_7(:,11);
D_7 = Mission_7(:,12);
P_7 = Mission_7(:,13);
E_7 = E_5b(end) -(Mission_7(1,14) - Mission_7(:,14)) ;  %%
Eperc_7 = (E_5b(end) -(Mission_7(1,14) - Mission_7(:,14)))/E_ini*100 ; %%


%% Deceleration

load('Mission_5c.mat')
t_5c = Mission_5c(:,1)+t_7(end); %%
x_5c = Mission_5c(:,2)+x_7(end);
h_5c = Mission_5c(:,3)+300;
V_5c = Mission_5c(:,4);
n_5c = Mission_5c(:,5);
gamma_5c = Mission_5c(:,6);
alpha_5c = Mission_5c(:,7);
epsilon_5c = Mission_5c(:,8);
phi_5c = Mission_5c(:,9);
L_5c = Mission_5c(:,10);
T_5c = Mission_5c(:,11);
D_5c = Mission_5c(:,12);
P_5c = Mission_5c(:,13);
E_5c = Mission_5c(:,14) - (E_ini-E_7(end)) ;  %%
Eperc_5c = Mission_5c(:,15) -(100 -Eperc_7(end)) ; %%
%% Descent

load('Mission_8.mat')
t_8 = Mission_8(:,1)+t_5c(end); %%
x_8 = Mission_8(:,2)+x_5c(end);
h_8 = Mission_8(:,3);
V_8 = Mission_8(:,4);
n_8 = Mission_8(:,5);
gamma_8 = Mission_8(:,6);
alpha_8 = Mission_8(:,7);
epsilon_8 = Mission_8(:,8);
phi_8 = Mission_8(:,9);
L_8 = Mission_8(:,10);
T_8 = Mission_8(:,11);
D_8 = Mission_8(:,12);
P_8 = Mission_8(:,13);
E_8 = E_5c(end) -Mission_8(:,14);  %%
Eperc_8 = Eperc_5c(end) -Mission_8(:,15) ; %%

%% H to V trans

load('Mission_9.mat')
t_9 = Mission_9(:,1)+t_8(end); %%
x_9 = Mission_9(:,2)+x_8(end);
h_9 = Mission_9(:,3);
V_9 = Mission_9(:,4);
n_9 = Mission_9(:,5);
gamma_9 = Mission_9(:,6);
alpha_9 = Mission_9(:,7);
epsilon_9 = Mission_9(:,8);
phi_9 = Mission_9(:,9);
L_9 = Mission_9(:,10);
T_9 = Mission_9(:,11);
D_9 = Mission_9(:,12);
P_9 = Mission_9(:,13);
E_9 = Mission_9(:,14) - (E_ini-E_8(end)) ;  %%
Eperc_9 = Mission_9(:,15) -(100 -Eperc_8(end)) ; %%

%% Hover 

load('Mission_10.mat')
t_10 = Mission_10(:,1)+t_9(end); %%
x_10 = Mission_10(:,2)+x_9(end);
h_10 = Mission_10(:,3)+50;
V_10 = Mission_10(:,4);
n_10 = Mission_10(:,5);
gamma_10 = Mission_10(:,6);
alpha_10 = Mission_10(:,7);
epsilon_10 = Mission_10(:,8);
phi_10 = Mission_10(:,9);
L_10 = Mission_10(:,10);
T_10 = Mission_10(:,11);
D_10 = Mission_10(:,12);
P_10 = Mission_10(:,13);
E_10 = Mission_10(:,14) - (E_ini-E_9(end)) ;  %%
Eperc_10 = Mission_10(:,15) -(100 -Eperc_9(end)) ; %%


%% Land

load('Mission_11.mat')
t_11 = Mission_11(:,1)+t_10(end); %%
x_11 = Mission_11(:,2)+x_10(end);
h_11 = Mission_11(:,3);
V_11 = Mission_11(:,4);
n_11 = Mission_11(:,5);
gamma_11 = Mission_11(:,6);
alpha_11 = Mission_11(:,7);
epsilon_11 = Mission_11(:,8);
phi_11 = Mission_11(:,9);
L_11 = Mission_11(:,10);
T_11 = Mission_11(:,11);
D_11 = Mission_11(:,12);
P_11 = Mission_11(:,13);
E_11 = Mission_11(:,14) - (E_ini-E_10(end)) ;  %%
Eperc_11 = Mission_11(:,15) -(100 -Eperc_10(end)) ; %%

%% Arrange vecs

t_vec = [t_1;t_2;t_3;t_4;t_5;t_6;t_5b;t_7;t_5c;t_8;t_9;t_10;t_11];
x_vec = [x_1;x_2;x_3;x_4;x_5;x_6;x_5b;x_7;x_5c;x_8;x_9;x_10;x_11];
h_vec = [h_1;h_2;h_3;h_4;h_5;h_6;h_5b;h_7;h_5c;h_8;h_9;h_10;h_11];
V_vec = [V_1;V_2;V_3;V_4;V_5;V_6;V_5b;V_7;V_5c;V_8;V_9;V_10;V_11];
n_vec = [n_1;n_2;n_3;n_4;n_5;n_6;n_5b;n_7;n_5c;n_8;n_9;n_10;n_11];
gamma_vec = [gamma_1;gamma_2;gamma_3;gamma_4;gamma_5;gamma_6;gamma_5b;gamma_7;gamma_5c;gamma_8;gamma_9;gamma_10;gamma_11];
alpha_vec = [alpha_1;alpha_2;alpha_3;alpha_4;alpha_5;alpha_6;alpha_5b;alpha_7;alpha_5c;alpha_8;alpha_9;alpha_10;alpha_11];
epsilon_vec = [epsilon_1;epsilon_2;epsilon_3;epsilon_4;epsilon_5;epsilon_6;epsilon_5b;epsilon_7;epsilon_5c;epsilon_8;epsilon_9;epsilon_10;epsilon_11];
phi_vec = [phi_1;phi_2;phi_3;phi_4;phi_5;phi_6;phi_5b;phi_7;phi_5c;phi_8;phi_9;phi_10;phi_11];
L_vec = [L_1;L_2;L_3;L_4;L_5;L_6;L_5b;L_7;L_5c;L_8;L_9;L_10;L_11];
T_vec = [T_1;T_2;T_3;T_4;T_5;T_6;T_5b;T_7;T_5c;T_8;T_9;T_10;T_11];
P_vec = [P_1;P_2;P_3;P_4;P_5;P_6;P_5b;P_7;P_5c;P_8;P_9;P_10;P_11];
D_vec = [D_1;D_2;D_3;D_4;D_5;D_6;D_5b;D_7;D_5c;D_8;D_9;D_10;D_11];
E_vec = [E_1;E_2;E_3;E_4;E_5;E_6;E_5b;E_7;E_5c;E_8;E_9;E_10;E_11];
Eperc_vec = [Eperc_1;Eperc_2;Eperc_3;Eperc_4;Eperc_5;Eperc_6;Eperc_5b;Eperc_7;Eperc_5c;Eperc_8;Eperc_9;Eperc_10;Eperc_11];



%% PLOTS!!!

figure(1)
plot(t_vec,x_vec,'b','LineWidth',1.2)
grid on
hold on
%yline(Pmax,'--b','LineWidth',1.2)
plot(t_1(end),x_1(end),'or','LineWidth',1.2)
plot(t_2(end),x_2(end),'or','LineWidth',1.2)
plot(t_3(end),x_3(end),'or','LineWidth',1.2)
plot(t_4(end),x_4(end),'or','LineWidth',1.2)
plot(t_5(end),x_5(end),'or','LineWidth',1.2)
plot(t_6(end),x_6(end),'or','LineWidth',1.2)
plot(t_5b(end),x_5b(end),'or','LineWidth',1.2)
plot(t_7(end),x_7(end),'or','LineWidth',1.2)
plot(t_5c(end),x_5c(end),'or','LineWidth',1.2)
plot(t_8(end),x_8(end),'or','LineWidth',1.2)
plot(t_9(end),x_9(end),'or','LineWidth',1.2)
plot(t_10(end),x_10(end),'or','LineWidth',1.2)
plot(t_11(end),x_11(end),'or','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('Horizontal distance [m]')
Title1 = strcat('Horizontal distance vs time. Full Mission');
title(Title1)

figure(2)
plot(t_vec,h_vec,'b','LineWidth',1.2)
grid on
hold on
%yline(Pmax,'--b','LineWidth',1.2)
plot(t_1(end),h_1(end),'or','LineWidth',1.2)
plot(t_2(end),h_2(end),'or','LineWidth',1.2)
plot(t_3(end),h_3(end),'or','LineWidth',1.2)
plot(t_4(end),h_4(end),'or','LineWidth',1.2)
plot(t_5(end),h_5(end),'or','LineWidth',1.2)
plot(t_6(end),h_6(end),'or','LineWidth',1.2)
plot(t_5b(end),h_5b(end),'or','LineWidth',1.2)
plot(t_7(end),h_7(end),'or','LineWidth',1.2)
plot(t_5c(end),h_5c(end),'or','LineWidth',1.2)
plot(t_8(end),h_8(end),'or','LineWidth',1.2)
plot(t_9(end),h_9(end),'or','LineWidth',1.2)
plot(t_10(end),h_10(end),'or','LineWidth',1.2)
plot(t_11(end),h_11(end),'or','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('Vertical distance [m]')
Title1 = strcat('Vertical distance vs time. Full Mission');
title(Title1)

figure(3)
plot(t_vec,V_vec,'b','LineWidth',1.2)
grid on
hold on
plot(t_1(end),V_1(end),'or','LineWidth',1.2)
plot(t_2(end),V_2(end),'or','LineWidth',1.2)
plot(t_3(end),V_3(end),'or','LineWidth',1.2)
plot(t_4(end),V_4(end),'or','LineWidth',1.2)
plot(t_5(end),V_5(end),'or','LineWidth',1.2)
plot(t_6(end),V_6(end),'or','LineWidth',1.2)
plot(t_5b(end),V_5b(end),'or','LineWidth',1.2)
plot(t_7(end),V_7(end),'or','LineWidth',1.2)
plot(t_5c(end),V_5c(end),'or','LineWidth',1.2)
plot(t_8(end),V_8(end),'or','LineWidth',1.2)
plot(t_9(end),V_9(end),'or','LineWidth',1.2)
plot(t_10(end),V_10(end),'or','LineWidth',1.2)
plot(t_11(end),V_11(end),'or','LineWidth',1.2)

%yline(Pmax,'--b','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('Total speed [m/s]')
Title1 = strcat('Aircraft speed vs time. Full Mission');
title(Title1)

figure(4)
plot(t_vec,n_vec,'b','LineWidth',1.2)
grid on
hold on
plot(t_1(end),n_1(end),'or','LineWidth',1.2)
plot(t_2(end),n_2(end),'or','LineWidth',1.2)
plot(t_3(end),n_3(end),'or','LineWidth',1.2)
plot(t_4(end),n_4(end),'or','LineWidth',1.2)
plot(t_5(end),n_5(end),'or','LineWidth',1.2)
plot(t_6(end),n_6(end),'or','LineWidth',1.2)
plot(t_5b(end),n_5b(end),'or','LineWidth',1.2)
plot(t_7(end),n_7(end),'or','LineWidth',1.2)
plot(t_5c(end),n_5c(end),'or','LineWidth',1.2)
plot(t_8(end),n_8(end),'or','LineWidth',1.2)
plot(t_9(end),n_9(end),'or','LineWidth',1.2)
plot(t_10(end),n_10(end),'or','LineWidth',1.2)
plot(t_11(end),n_11(end),'or','LineWidth',1.2)
yline(rps_max,'--k','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('Engine revolutions [rps]')
Title1 = strcat('Engine revolutions vs time. Full Mission');
title(Title1)


figure(5)
plot(t_vec,gamma_vec*180/pi,'b','LineWidth',1.2)
grid on
hold on
plot(t_1(end),gamma_1(end)*180/pi,'or','LineWidth',1.2)
plot(t_2(end),gamma_2(end)*180/pi,'or','LineWidth',1.2)
plot(t_3(end),gamma_3(end)*180/pi,'or','LineWidth',1.2)
plot(t_4(end),gamma_4(end)*180/pi,'or','LineWidth',1.2)
plot(t_5(end),gamma_5(end)*180/pi,'or','LineWidth',1.2)
plot(t_6(end),gamma_6(end)*180/pi,'or','LineWidth',1.2)
plot(t_5b(end),gamma_5b(end)*180/pi,'or','LineWidth',1.2)
plot(t_7(end),gamma_7(end)*180/pi,'or','LineWidth',1.2)
plot(t_5c(end),gamma_5c(end)*180/pi,'or','LineWidth',1.2)
plot(t_8(end),gamma_8(end)*180/pi,'or','LineWidth',1.2)
plot(t_9(end),gamma_9(end)*180/pi,'or','LineWidth',1.2)
plot(t_10(end),gamma_10(end)*180/pi,'or','LineWidth',1.2)
plot(t_11(end),gamma_11(end)*180/pi,'or','LineWidth',1.2)
%yline(Pmax,'--b','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('\gamma [deg]')
Title1 = strcat('Trajectory angle vs time. Full Mission');
title(Title1)

figure(6)
plot(t_vec,alpha_vec*180/pi,'b','LineWidth',1.2)
grid on
hold on
plot(t_1(end),alpha_1(end)*180/pi,'or','LineWidth',1.2)
plot(t_2(end),alpha_2(end)*180/pi,'or','LineWidth',1.2)
plot(t_3(end),alpha_3(end)*180/pi,'or','LineWidth',1.2)
plot(t_4(end),alpha_4(end)*180/pi,'or','LineWidth',1.2)
plot(t_5(end),alpha_5(end)*180/pi,'or','LineWidth',1.2)
plot(t_6(end),alpha_6(end)*180/pi,'or','LineWidth',1.2)
plot(t_5b(end),alpha_5b(end)*180/pi,'or','LineWidth',1.2)
plot(t_7(end),alpha_7(end)*180/pi,'or','LineWidth',1.2)
plot(t_5c(end),alpha_5c(end)*180/pi,'or','LineWidth',1.2)
plot(t_8(end),alpha_8(end)*180/pi,'or','LineWidth',1.2)
plot(t_9(end),alpha_9(end)*180/pi,'or','LineWidth',1.2)
plot(t_10(end),alpha_10(end)*180/pi,'or','LineWidth',1.2)
plot(t_11(end),alpha_11(end)*180/pi,'or','LineWidth',1.2)
%yline(Pmax,'--b','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('\alpha [deg]')
Title1 = strcat('Angle of attack vs time. Full Mission');
title(Title1)

figure(7)
plot(t_vec,epsilon_vec*180/pi,'b','LineWidth',1.2)
grid on
hold on
plot(t_1(end),epsilon_1(end)*180/pi,'or','LineWidth',1.2)
plot(t_2(end),epsilon_2(end)*180/pi,'or','LineWidth',1.2)
plot(t_3(end),epsilon_3(end)*180/pi,'or','LineWidth',1.2)
plot(t_4(end),epsilon_4(end)*180/pi,'or','LineWidth',1.2)
plot(t_5(end),epsilon_5(end)*180/pi,'or','LineWidth',1.2)
plot(t_6(end),epsilon_6(end)*180/pi,'or','LineWidth',1.2)
plot(t_5b(end),epsilon_5b(end)*180/pi,'or','LineWidth',1.2)
plot(t_7(end),epsilon_7(end)*180/pi,'or','LineWidth',1.2)
plot(t_5c(end),epsilon_5c(end)*180/pi,'or','LineWidth',1.2)
plot(t_8(end),epsilon_8(end)*180/pi,'or','LineWidth',1.2)
plot(t_9(end),epsilon_9(end)*180/pi,'or','LineWidth',1.2)
plot(t_10(end),epsilon_10(end)*180/pi,'or','LineWidth',1.2)
plot(t_11(end),epsilon_11(end)*180/pi,'or','LineWidth',1.2)
%yline(Pmax,'--b','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('\epsilon [deg]')
Title1 = strcat('Engine tilt angle vs time. Full Mission');
title(Title1)

figure(8)
plot(t_vec,phi_vec*180/pi,'b','LineWidth',1.2)
grid on
hold on
plot(t_1(end),phi_1(end)*180/pi,'or','LineWidth',1.2)
plot(t_2(end),phi_2(end)*180/pi,'or','LineWidth',1.2)
plot(t_3(end),phi_3(end)*180/pi,'or','LineWidth',1.2)
plot(t_4(end),phi_4(end)*180/pi,'or','LineWidth',1.2)
plot(t_5(end),phi_5(end)*180/pi,'or','LineWidth',1.2)
plot(t_6(end),phi_6(end)*180/pi,'or','LineWidth',1.2)
plot(t_5b(end),phi_5b(end)*180/pi,'or','LineWidth',1.2)
plot(t_7(end),phi_7(end)*180/pi,'or','LineWidth',1.2)
plot(t_5c(end),phi_5c(end)*180/pi,'or','LineWidth',1.2)
plot(t_8(end),phi_8(end)*180/pi,'or','LineWidth',1.2)
plot(t_9(end),phi_9(end)*180/pi,'or','LineWidth',1.2)
plot(t_10(end),phi_10(end)*180/pi,'or','LineWidth',1.2)
plot(t_11(end),phi_11(end)*180/pi,'or','LineWidth',1.2)
%yline(Pmax,'--b','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('\phi [deg]')
Title1 = strcat('Total inclination angle vs time. Full Mission');
title(Title1)

figure(9)
plot(t_vec,L_vec,'b','LineWidth',1.2)
grid on
hold on
plot(t_1(end),L_1(end),'or','LineWidth',1.2)
plot(t_2(end),L_2(end),'or','LineWidth',1.2)
plot(t_3(end),L_3(end),'or','LineWidth',1.2)
plot(t_4(end),L_4(end),'or','LineWidth',1.2)
plot(t_5(end),L_5(end),'or','LineWidth',1.2)
plot(t_6(end),L_6(end),'or','LineWidth',1.2)
plot(t_5b(end),L_5b(end),'or','LineWidth',1.2)
plot(t_7(end),L_7(end),'or','LineWidth',1.2)
plot(t_5c(end),L_5c(end),'or','LineWidth',1.2)
plot(t_8(end),L_8(end),'or','LineWidth',1.2)
plot(t_9(end),L_9(end),'or','LineWidth',1.2)
plot(t_10(end),L_10(end),'or','LineWidth',1.2)
plot(t_11(end),L_11(end),'or','LineWidth',1.2)
%yline(Pmax,'--b','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('L [N]')
Title1 = strcat('Aerodynamic Lift force vs time. Full Mission');
title(Title1)

figure(10)
plot(t_vec,T_vec,'b','LineWidth',1.2)
grid on
hold on
plot(t_1(end),T_1(end),'or','LineWidth',1.2)
plot(t_2(end),T_2(end),'or','LineWidth',1.2)
plot(t_3(end),T_3(end),'or','LineWidth',1.2)
plot(t_4(end),T_4(end),'or','LineWidth',1.2)
plot(t_5(end),T_5(end),'or','LineWidth',1.2)
plot(t_6(end),T_6(end),'or','LineWidth',1.2)
plot(t_5b(end),T_5b(end),'or','LineWidth',1.2)
plot(t_7(end),T_7(end),'or','LineWidth',1.2)
plot(t_5c(end),T_5c(end),'or','LineWidth',1.2)
plot(t_8(end),T_8(end),'or','LineWidth',1.2)
plot(t_9(end),T_9(end),'or','LineWidth',1.2)
plot(t_10(end),T_10(end),'or','LineWidth',1.2)
plot(t_11(end),T_11(end),'or','LineWidth',1.2)
yline(Tmax,'--b','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('T [N]')
Title1 = strcat('Thrust force vs time. Full Mission');
title(Title1)

figure(11)
plot(t_vec,D_vec,'b','LineWidth',1.2)
grid on
hold on
plot(t_1(end),D_1(end),'or','LineWidth',1.2)
plot(t_2(end),D_2(end),'or','LineWidth',1.2)
plot(t_3(end),D_3(end),'or','LineWidth',1.2)
plot(t_4(end),D_4(end),'or','LineWidth',1.2)
plot(t_5(end),D_5(end),'or','LineWidth',1.2)
plot(t_6(end),D_6(end),'or','LineWidth',1.2)
plot(t_5b(end),D_5b(end),'or','LineWidth',1.2)
plot(t_7(end),D_7(end),'or','LineWidth',1.2)
plot(t_5c(end),D_5c(end),'or','LineWidth',1.2)
plot(t_8(end),D_8(end),'or','LineWidth',1.2)
plot(t_9(end),D_9(end),'or','LineWidth',1.2)
plot(t_10(end),D_10(end),'or','LineWidth',1.2)
plot(t_11(end),D_11(end),'or','LineWidth',1.2)
%yline(Pmax,'--b','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('D [N]')
Title1 = strcat('Aerodynamic Drag force vs time. Full Mission');
title(Title1)

figure(12)
plot(t_vec,P_vec,'b','LineWidth',1.2)
grid on
hold on
plot(t_1(end),P_1(end),'or','LineWidth',1.2)
plot(t_2(end),P_2(end),'or','LineWidth',1.2)
plot(t_3(end),P_3(end),'or','LineWidth',1.2)
plot(t_4(end),P_4(end),'or','LineWidth',1.2)
plot(t_5(end),P_5(end),'or','LineWidth',1.2)
plot(t_6(end),P_6(end),'or','LineWidth',1.2)
plot(t_5b(end),P_5b(end),'or','LineWidth',1.2)
plot(t_7(end),P_7(end),'or','LineWidth',1.2)
plot(t_5c(end),P_5c(end),'or','LineWidth',1.2)
plot(t_8(end),P_8(end),'or','LineWidth',1.2)
plot(t_9(end),P_9(end),'or','LineWidth',1.2)
plot(t_10(end),P_10(end),'or','LineWidth',1.2)
plot(t_11(end),P_11(end),'or','LineWidth',1.2)
%yline(Pmax,'--b','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('P [W]')
Title1 = strcat('Power consumption vs time. Full Mission');
title(Title1)

figure(13)
plot(t_vec,E_vec,'b','LineWidth',1.2)
grid on
hold on
plot(t_1(end),E_1(end),'or','LineWidth',1.2)
plot(t_2(end),E_2(end),'or','LineWidth',1.2)
plot(t_3(end),E_3(end),'or','LineWidth',1.2)
plot(t_4(end),E_4(end),'or','LineWidth',1.2)
plot(t_5(end),E_5(end),'or','LineWidth',1.2)
plot(t_6(end),E_6(end),'or','LineWidth',1.2)
plot(t_5b(end),E_5b(end),'or','LineWidth',1.2)
plot(t_7(end),E_7(end),'or','LineWidth',1.2)
plot(t_5c(end),E_5c(end),'or','LineWidth',1.2)
plot(t_8(end),E_8(end),'or','LineWidth',1.2)
plot(t_9(end),E_9(end),'or','LineWidth',1.2)
plot(t_10(end),E_10(end),'or','LineWidth',1.2)
plot(t_11(end),E_11(end),'or','LineWidth',1.2)
yline(e0*mbatt*0.2,'--k','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('E [J]')
Title1 = strcat('Energy consumption vs time. Full Mission');
title(Title1)

figure(14)
plot(t_vec,Eperc_vec,'b','LineWidth',1.2)
grid on
hold on
plot(t_1(end),Eperc_1(end),'or','LineWidth',1.2)
plot(t_2(end),Eperc_2(end),'or','LineWidth',1.2)
plot(t_3(end),Eperc_3(end),'or','LineWidth',1.2)
plot(t_4(end),Eperc_4(end),'or','LineWidth',1.2)
plot(t_5(end),Eperc_5(end),'or','LineWidth',1.2)
plot(t_6(end),Eperc_6(end),'or','LineWidth',1.2)
plot(t_5b(end),Eperc_5b(end),'or','LineWidth',1.2)
plot(t_7(end),Eperc_7(end),'or','LineWidth',1.2)
plot(t_5c(end),Eperc_5c(end),'or','LineWidth',1.2)
plot(t_8(end),Eperc_8(end),'or','LineWidth',1.2)
plot(t_9(end),Eperc_9(end),'or','LineWidth',1.2)
plot(t_10(end),Eperc_10(end),'or','LineWidth',1.2)
plot(t_11(end),Eperc_11(end),'or','LineWidth',1.2)
yline(20,'--k','LineWidth',1.2)
xlabel('Mission time [s] ')
ylabel('E [%]')
Title1 = strcat('Energy consumption (percentage) vs time. Full Mission');
title(Title1)