function transition_v6
%Evolutions epsilon(V) and alpha(V) are considered.

%Initial data:
m        = 16.6;            %[kg]
g        = 9.81;            %[m/s^2]
S        = 0.5008;          %[m^2]
rho      = 1.2133;          %[kg/m^3]
alpha_mp = 16.4*pi/180;     %[rad]
alpha_mi = -5*pi/180;       %[rad]
N_eng    = 2;               %[-]
D_prop   = 0.8128;          %[m]
eta_m    = 0.73;            %[-]
n_max    = 78.125;          %[rps]
tau_tr   = 8;              %[s]
alpha_t  = 12.9217*pi/180;  %[rad]
V_f      = 20;              %[m/s]

%1) Heuristic transition:
params.m        = m;
params.g        = g;
params.S        = S;
params.rho      = rho;
params.N_eng    = N_eng;
params.D_prop   = D_prop;
params.eta_m    = eta_m;
params.n_max    = n_max;
params.alpha_mi = alpha_mi;
params.alpha_t  = alpha_t;
params.V_f      = V_f;

[CL_mp,CD_mp] = aerodyn(alpha_mp);
V_mp          = sqrt(2*m*g/(rho*S*(CL_mp+CD_mp*tan(alpha_mp))));

%Initial state:
[C_T0,C_P0] = propmodel(0,0);
n_hov       = sqrt(m*g/(N_eng*rho*D_prop^4*C_T0));
P_0         = N_eng*rho*D_prop^5*n_hov^3*C_P0;

y_p0    = [0;P_0/(1000*eta_m);0]; %Resulting energy is in kJ!!
y_0     = [0;0;n_hov];

ode_ops = odeset('AbsTol',1e-11,'RelTol',1e-11,'Events',@(t,y) event_nmax(t,y,n_max),...
          'Mass',[1 0 0; 0 1 0; 0 0 0],'MStateDependence','none','MassSingular','yes',...
          'InitialSlope',y_p0);
trans_1 = ode15s(@(t,y) dynfunDAE_heuris(t,y,params),[0;tau_tr],y_0,ode_ops);

t_vec    = trans_1.x;
V_vec    = trans_1.y(1,:);
E_vec    = trans_1.y(2,:);
n_vec    = trans_1.y(3,:);
alp_vec  = alpha_mi*(1-exp(-4*t_vec)) + (alpha_t-alpha_mi)*...
    (exp(min([V_vec; V_f*ones(size(V_vec))],[],1)/4)./exp(V_f*ones(size(V_vec))/4));
epsi_vec = 0.1*pi/2*exp(-4*t_vec) + ...
    0.9*pi/2*(1-exp(min([V_vec; V_f*ones(size(V_vec))],[],1)/4)./exp(V_f*ones(size(V_vec))/4));
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
P_vec           = N_eng*rho*D_prop^5*n_vec.^3.*C_Pvec;

figure(5), hold on, grid on,
plot(t_vec,T_vec,'k'),plot(t_vec,L_vec,'b'), plot(t_vec,D_vec,'r')
xlabel('t [s]'), ylabel('T, L, and D [N]'), title('Forces')

figure(6), hold on, grid on, 
plot(V_vec,180/pi*epsi_vec,'k')
xlabel('V [m/s]'), ylabel('epsilon [deg]'), title('Tilt angle trajectory')

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
N_eng    = params.N_eng;
D_prop   = params.D_prop;
eta_m    = params.eta_m;
alpha_mi = params.alpha_mi;
alpha_t  = params.alpha_t;
V_f      = params.V_f;

%Variables:
V = y(1);
n = y(3);

%Computations:
J    = V/(n*D_prop);
if V<V_f
%     alp  = alpha_mi + (alpha_t-alpha_mi)*(3*(V/V_f)^2-2*(V/V_f)^3);
    alp  = alpha_mi*(1-exp(-4*t)) + (alpha_t-alpha_mi)*(exp(V/4)/exp(V_f/4));
    epsi = 0.1*pi/2*exp(-4*t) + 0.9*pi/2*(1-(exp(V/4)/exp(V_f/4)));
else
    alp  = alpha_t;
    epsi = 0;
end
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



