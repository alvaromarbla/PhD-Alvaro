function transition_v3

%Initial data:
m        = 16.6;            %[kg]
g        = 9.81;            %[m/s^2]
S        = 0.5008;          %[m^2]
rho      = 1.2133;          %[kg/m^3]
alpha_mp = 16.4*pi/180;     %[rad]
alpha_t  = 12.9217*pi/180;  %[rad]
alpha_mi = -0*pi/180;       %[rad]
N_eng    = 2;               %[-]
D_prop   = 0.8128;          %[m]
eta_m    = 0.88;            %[-]
n_max    = 78.125;          %[rps]
tau_tr   = 12;              %[s]
C_cru    = 16.7;            %[J/m]


%1) Heuristic transition:
[CL_mp,CD_mp] = aerodyn(alpha_mp);
V_mp          = sqrt(2*m*g/(rho*S*CL_mp));
[C_T0,C_P0]   = propmodel(0,0);
n_hov         = sqrt(m*g/(N_eng*rho*D_prop^4*C_T0));
P_0           = N_eng*rho*D_prop^5*n_hov^3*C_P0;

params.m        = m;
params.g        = g;
params.S        = S;
params.rho      = rho;
params.alpha_mp = alpha_mp;
params.alpha_t  = alpha_t;
params.alpha_mi = alpha_mi;
params.N_eng    = N_eng;
params.D_prop   = D_prop;
params.eta_m    = eta_m;
params.n_max    = n_max;
params.tau_tr   = tau_tr;
params.C_cru    = C_cru;

y_p0    = [0;P_0/(1000*eta_m);0]; %Resulting energy is in kJ!!
y_0     = [0;0;n_hov];

ode_ops = odeset('AbsTol',1e-11,'RelTol',1e-11,'Events',@(t,y) event_nmax(t,y,n_max),...
          'Mass',[1 0 0; 0 1 0; 0 0 0],'MStateDependence','none','MassSingular','yes',...
          'InitialSlope',y_p0);
trans_1 = ode15s(@(t,y) dynfunDAE_heuris(t,y,params),[0;tau_tr],y_0,ode_ops);

t_vec    = trans_1.x;

%Alpha(t): Opciones:
% alp_vec  = alpha_mi + (alpha_mp-alpha_mi)/tau_tr*t_vec;
% alp_vec  = alpha_mi + (alpha_mp/2-alpha_mi)/tau_tr*t_vec;
% alp_vec  = alpha_mi - 2*(alpha_t-alpha_mi)*((t_vec/tau_tr).^3-1.5*(t_vec/tau_tr).^2);
% alp_vec = alpha_t + (alpha_mi-alpha_t)*(1-t_vec/tau_tr).^3;
alp_vec  = alpha_t + (alpha_mi-alpha_t)*exp(-5*t_vec/tau_tr);

V_vec    = trans_1.y(1,:);
E_vec    = trans_1.y(2,:);
n_vec    = trans_1.y(3,:);

%Epsilon(t): Opciones:
% epsi_vec = pi/2*(1-(t_vec/tau_tr).^3);
% epsi_vec = pi/2*(1-min([V_vec/20; ones(size(V_vec))],[],1).^3);
% epsi_vec = pi/2*(1-(t_vec/tau_tr));
% epsi_vec = pi/2*exp(-5*t_vec/tau_tr);
% epsi_vec = pi/2 + pi*((t_vec/tau_tr).^3-1.5*(t_vec/tau_tr).^2);
epsi_vec = pi/2*(1-t_vec/tau_tr).^2;

phi_vec  = epsi_vec + alp_vec;

figure(1), hold on, grid on, plot(t_vec,V_vec)
plot(trans_1.x,V_mp*ones(size(trans_1.x)),'-.')
xlabel('t [s]'), ylabel('V [m/s]'), title('Airspeed')

figure(2), hold on, grid on, plot(t_vec,E_vec)
xlabel('t [s]'), ylabel('\Delta E [kJ]'), title('Energy')

figure(3), hold on, grid on, plot(t_vec,n_vec)
xlabel('t [s]'), ylabel('n [rps]'), title('Engine running speed')

figure(4), hold on, grid on, plot(t_vec,180/pi*phi_vec)
plot(t_vec,180/pi*epsi_vec,'k'),plot(t_vec,180/pi*alp_vec,'r')
xlabel('t [s]'), ylabel('phi [deg]'), title('Thrust angle of attack')

[CL_vec,CD_vec] = aerodyn(alp_vec);
L_vec = 0.5*rho*V_vec.^2*S.*CL_vec;
D_vec = 0.5*rho*V_vec.^2*S.*CD_vec;

J_vec           = V_vec./(D_prop*n_vec);
[C_Tvec,C_Pvec] = propmodel(J_vec,phi_vec);
T_vec           = N_eng*rho*D_prop^4*n_vec.^2.*C_Tvec;
% P_vec           = N_eng*rho*D_prop^5*n_vec.^3.*C_Pvec;

figure(5), hold on, grid on,
plot(t_vec,T_vec,'k'),plot(t_vec,L_vec,'b'), plot(t_vec,D_vec,'r')
xlabel('t [s]'), ylabel('T, L, and D [N]'), title('Forces')


%Sensitivity analysis:
tau_tr          = 0.2;
params.epsi_s   = tau_tr;

trans_2 = ode15s(@(t,y) dynfunDAE_heuris(t,y,params),[0;(pi/2-alpha_mp)/tau_tr],y_0,ode_ops);

t_vec   = trans_2.x;
phi_vec = (pi/2 - tau_tr*t_vec);
V_vec   = trans_2.y(1,:);
E_vec   = trans_2.y(2,:);
n_vec   = trans_2.y(3,:);

L_vec = 0.5*rho*V_vec.^2*S*CL_mp;
D_vec = 0.5*rho*V_vec.^2*S*CD_mp;

J_vec           = V_vec./(D_prop*n_vec);
[C_Tvec,C_Pvec] = propmodel(J_vec,phi_vec);
T_vec           = N_eng*rho*D_prop^4*n_vec.^2.*C_Tvec;
% P_vec           = N_eng*rho*D_prop^5*n_vec.^3.*C_Pvec;

figure(1), hold on, grid on, plot(t_vec,V_vec,'r')
figure(2), hold on, grid on, plot(t_vec,E_vec,'r')
figure(3), hold on, grid on, plot(t_vec,n_vec,'r')
figure(4), hold on, grid on, plot(t_vec,180/pi*phi_vec,'r')
figure(5), plot(t_vec,T_vec,'k--'),plot(t_vec,L_vec,'b--'), plot(t_vec,D_vec,'r--')

%Sensitivity analysis:
tau_tr          = 0.4;
params.epsi_s   = tau_tr;

trans_3 = ode15s(@(t,y) dynfunDAE_heuris(t,y,params),[0;(pi/2-alpha_mp)/tau_tr],y_0,ode_ops);

t_vec   = trans_3.x;
phi_vec = (pi/2 - tau_tr*t_vec);
V_vec   = trans_3.y(1,:);
E_vec   = trans_3.y(2,:);
n_vec   = trans_3.y(3,:);

L_vec = 0.5*rho*V_vec.^2*S*CL_mp;
D_vec = 0.5*rho*V_vec.^2*S*CD_mp;

J_vec           = V_vec./(D_prop*n_vec);
[C_Tvec,C_Pvec] = propmodel(J_vec,phi_vec);
T_vec           = N_eng*rho*D_prop^4*n_vec.^2.*C_Tvec;
% P_vec           = N_eng*rho*D_prop^5*n_vec.^3.*C_Pvec;

figure(1), hold on, grid on, plot(t_vec,V_vec,'k')
figure(2), hold on, grid on, plot(t_vec,E_vec,'k')
figure(3), hold on, grid on, plot(t_vec,n_vec,'k')
figure(4), hold on, grid on, plot(t_vec,180/pi*phi_vec,'k')
figure(5), plot(t_vec,T_vec,'k-.'),plot(t_vec,L_vec,'b-.'), plot(t_vec,D_vec,'r-.')


disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = dynfunDAE_heuris(t,y,params)
%This function implements the a/c dynamics with y=[V;E;n]. n is treated as
%a dynamic variable, but its evolution comes from the algebraic equation at
%each time instant.

%Parameters:
m        = params.m;
g        = params.g;
S        = params.S;
rho      = params.rho;
alpha_mp = params.alpha_mp;
alpha_t  = params.alpha_t;
alpha_mi = params.alpha_mi;
N_eng    = params.N_eng;
D_prop   = params.D_prop;
eta_m    = params.eta_m;
tau_tr   = params.tau_tr;

%Variables:
V = y(1);
n = y(3);

%Computations:
J    = V/(n*D_prop);

%Options of alp:
% alp  = alpha_mi + (alpha_mp-alpha_mi)/tau_tr*t;
% alp  = alpha_mi + (alpha_mp-alpha_mi)*(t/tau_tr)^2;
% alp  = alpha_mi + (alpha_mp/2-alpha_mi)/tau_tr*t;
% alp  = 0;
% alp  = alpha_mi - 2*(alpha_t-alpha_mi)*((t/tau_tr).^3-1.5*(t/tau_tr).^2);
% alp = alpha_t + (alpha_mi-alpha_t)*(1-t/tau_tr).^3;
alp  = alpha_t + (alpha_mi-alpha_t)*exp(-5*t/tau_tr);

%Options of epsi:
% epsi = pi/2-(pi/2+alpha_mp)*t/tau_tr;
% epsi = pi/2*(1-(t/tau_tr)^3);
% epsi = pi/2*(1-min([V/20;1])^3); 
% epsi = pi/2*(1-(t/tau_tr));
% epsi = pi/2*exp(-5*t/tau_tr);
% epsi = pi/2 + pi*((t/tau_tr).^3-1.5*(t/tau_tr).^2);
epsi = pi/2*(1-t/tau_tr)^2;

phi  = epsi + alp;

[C_T,C_P] = propmodel(J,phi);
T         = N_eng*rho*D_prop^4*n^2*C_T;
P         = N_eng*rho*D_prop^5*n^3*C_P;

[C_L,C_D] = aerodyn(alp);
L         = 0.5*rho*V^2*S*C_L;
D         = 0.5*rho*V^2*S*C_D;

f    = zeros(3,1);
f(1) = (T*cos(phi)-D)/m;
f(2) = P/(1000*eta_m);              %Resulting energy is in kJ!!
f(3) = (L + T*sin(phi))/(m*g)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value,isterminal,direction] = event_nmax(t,y,n_max)

value      = y(3)-n_max;
isterminal = 1;
direction  = 0;



