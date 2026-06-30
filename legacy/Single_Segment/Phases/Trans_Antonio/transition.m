function transition

%Initial data:
m        = 16.6;        %[kg]
g        = 9.81;        %[m/s^2]
S        = 0.5008;      %[m^2]
rho      = 1.2133;      %[kg/m^3]
alpha_mp = 2.6*pi/180; %[rad]
N_eng    = 2;           %[-]
D_prop   = 0.8128;      %[m]
eta_m    = 0.88;        %[-]
n_max    = 78.125;      %[rps]
epsi_s   = 0.3;         %[rad/s] % this was 0.1 before
C_cru    = 16.7;        %[J/m]
SF       = 1.2;

%1) Heuristic transition:
[CL_mp,CD_mp] = aerodyn(alpha_mp);
CL_Max        = 2.045;
CL_Max_ope    = CL_Max/SF^2;
V_mp          = sqrt(2*m*g/(rho*S*CL_Max_ope));
[C_T0,C_P0]   = propmodel(0,0);
n_hov         = sqrt(m*g/(N_eng*rho*D_prop^4*C_T0));
P_0           = N_eng*rho*D_prop^5*n_hov^3*C_P0;

params.m        = m;
params.g        = g;
params.S        = S;
params.rho      = rho;
params.alpha_mp = alpha_mp;
params.N_eng    = N_eng;
params.D_prop   = D_prop;
params.eta_m    = eta_m;
params.n_max    = n_max;
params.epsi_s   = epsi_s;
params.C_cru    = C_cru;

y_p0    = [0;P_0/(1000*eta_m);0]; %Resulting energy is in kJ!!
y_0     = [0;0;n_hov];

ode_ops = odeset('AbsTol',1e-11,'RelTol',1e-11,'Events',@(t,y) event_nmax(t,y,V_mp),...
          'Mass',[1 0 0; 0 1 0; 0 0 0],'MStateDependence','none','MassSingular','yes',...
          'InitialSlope',y_p0);
trans_1 = ode15s(@(t,y) dynfunDAE_almax(t,y,params),[0;(pi/2-alpha_mp)/epsi_s],y_0,ode_ops);

t_vec   = trans_1.x;
phi_vec = (pi/2 - epsi_s*t_vec);
V_vec   = trans_1.y(1,:);
E_vec   = trans_1.y(2,:);
n_vec   = trans_1.y(3,:);

figure(1), hold on, grid on, plot(t_vec,V_vec)
plot(trans_1.x,V_mp*ones(size(trans_1.x)),'-.')
xlabel('t [s]'), ylabel('V [m/s]'), title('Airspeed')

figure(2), hold on, grid on, plot(t_vec,E_vec)
xlabel('t [s]'), ylabel('\Delta E [kJ]'), title('Energy')

figure(3), hold on, grid on, plot(t_vec,n_vec)
xlabel('t [s]'), ylabel('n [rps]'), title('Engine running speed')

figure(4), hold on, grid on, plot(t_vec,180/pi*phi_vec)
xlabel('t [s]'), ylabel('phi [deg]'), title('Thrust angle of attack')

L_vec = 0.5*rho*V_vec.^2*S*CL_mp;
D_vec = 0.5*rho*V_vec.^2*S*CD_mp;

J_vec           = V_vec./(D_prop*n_vec);
[C_Tvec,C_Pvec] = propmodel(J_vec,phi_vec);
T_vec           = N_eng*rho*D_prop^4*n_vec.^2.*C_Tvec;
% P_vec           = N_eng*rho*D_prop^5*n_vec.^3.*C_Pvec;

figure(5), hold on, grid on,
plot(t_vec,T_vec,'k'),plot(t_vec,L_vec,'b'), plot(t_vec,D_vec,'r')
xlabel('t [s]'), ylabel('T, L, and D [N]'), title('Forces')


%Sensitivity analysis:
epsi_s          = 0.3; % this was 0.2 before
params.epsi_s   = epsi_s;

trans_2 = ode15s(@(t,y) dynfunDAE_almax(t,y,params),[0;(pi/2-alpha_mp)/epsi_s],y_0,ode_ops);

t_vec   = trans_2.x;
phi_vec = (pi/2 - epsi_s*t_vec);
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
epsi_s          = 0.3; % this was 0.4 before 
params.epsi_s   = epsi_s;

trans_3 = ode15s(@(t,y) dynfunDAE_almax(t,y,params),[0;(pi/2-alpha_mp)/epsi_s],y_0,ode_ops);

t_vec   = trans_3.x;
phi_vec = (pi/2 - epsi_s*t_vec);
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
function f = dynfunDAE_almax(t,y,params)
%This function implements the a/c dynamics with y=[V;E;n]. n is treated as
%a dynamic variable, but its evolution comes from the algebraic equation at
%each time instant.

%Parameters:
m        = params.m;
g        = params.g;
S        = params.S;
rho      = params.rho;
alpha_mp = params.alpha_mp;
N_eng    = params.N_eng;
D_prop   = params.D_prop;
eta_m    = params.eta_m;
epsi_s   = params.epsi_s;

%Variables:
V = y(1);
n = y(3);

%Computations:
J   = V/(n*D_prop);
phi = pi/2 - epsi_s*t;

[C_T,C_P] = propmodel(J,phi);
T         = N_eng*rho*D_prop^4*n^2*C_T;
P         = N_eng*rho*D_prop^5*n^3*C_P;

[C_L,C_D] = aerodyn(alpha_mp);
L         = 0.5*rho*V^2*S*C_L;
D         = 0.5*rho*V^2*S*C_D;

f    = zeros(3,1);
f(1) = (T*cos(phi)-D)/m;
f(2) = P/(1000*eta_m);              %Resulting energy is in kJ!!
% f(3) = (L + T*sin(phi))/(m*g)-1;
f(3) = (L + T*sin(phi))/(m*g)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value,isterminal,direction] = event_nmax(t,y,V_mp)
%n_max
value      = y(1) - V_mp; % y(3)-n_max;
isterminal = 1;
direction  = 0;



