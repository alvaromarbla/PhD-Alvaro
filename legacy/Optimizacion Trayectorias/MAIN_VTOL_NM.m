% MAIN.m -- Goddard Rocket
%
% This script runs a trajectory optimization to find the optimal thrust
% trajectory for the EMERGENTIA EVTOL RPAS to reach a desired altitude. ¡
%
% Dynamics come with a speed dependant model for drag and a wind tunnel
% model for thrust and power.
%

clc; clear;
addpath ../../

%% Variables of the problem: 

D     = 0.8128; % Propeller Diameter [m]
rho   = 1.2133; % Air density [kg/m^3]
N_eng = 2; % Number of engines

deltaH = 50; %[m]
%% Define maximum rps

RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;

%% Loads and weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
S    = 0.5008; % Reference surface [m^2]

%% Electric
tau = 0.2;
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency
mbatt = 3; %[kg]
e0    = 720e3; %[J/kg]
E_batt     = e0*mbatt;  %[J] Total energy of the battery packs
E_empty    = E_batt*tau;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Problem Bounds                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

h0 = 0;  %Vehicle starts on the ground
v0 = 0;  %Vehicle starts stationary
E0 = E_batt;  %Vehicle starts full of battery

hF = deltaH; % Vehicle needs to reach a height
vF = 0;  %Vehicle needs to be stationary at the end
EF = E_batt*tau; % Minimum energy we can ensure

hLow = 0;   %Cannot go through the earth
hUpp = deltaH*1.05;  %Maybe some overoscillation

vLow = 0; %Just look at the trajectory as it goes up
vUpp = 15;  % Go as fast as you can

nHover = 45.7558; % Result from Hovering
uLow = nHover/2; %%40
uUpp = rps_max; %Maximum rps output

% PLow = 0;
% PUpp = Pmax; %Maximum Power output

P.bounds.initialTime.low = 0;
P.bounds.initialTime.upp = 0;

P.bounds.finalTime.low = 0;
P.bounds.finalTime.upp = 60;

P.bounds.state.low = [hLow;vLow;EF];
P.bounds.state.upp = [hUpp;vUpp;E0];

P.bounds.initialState.low = [h0;v0;E0];
P.bounds.initialState.upp = [h0;v0;E0];

P.bounds.finalState.low = [hF;vF;EF];
P.bounds.finalState.upp = [hF;vF;E0];

P.bounds.control.low = uLow;
P.bounds.control.upp = uUpp;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Initial Guess                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
hGuess = hF;   %(m) guess at the maximum height reached
nHover = 45.7558; % Result from Hovering
nIni   = 50; % Result from Optimal ode45
P.guess.time    = [0, 20];  %(s)
P.guess.state   = [ [h0;v0;E0],  [hGuess;vF;E0*0.94] ];
P.guess.control = [nIni nHover];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Objective and Dynamic functions                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Dynamics function:
P.func.dynamics = @(t,x,u)( VTOLDynamics_NM(x,u) );

% Objective function:
P.func.bndObj = @(t0,x0,tF,xF)( xF(3));  %Minimize final energy /(tF-t0) + tF

%P.func.pathObj = @(t,x,u)( (u-nIni).^2);  %Minimize control output

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  Options and Method selection                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%method = 'trapezoid';
method = 'rungeKutta';
%method = 'chebyshev';
%method = 'hermiteSimpson';

switch method
    
    case 'trapezoid'
        
        P.options(1).method = 'trapezoid';
        P.options(1).defaultAccuracy = 'low';
        
        P.options(2).method = 'trapezoid';
        P.options(2).defaultAccuracy = 'medium';
        P.options(2).nlpOpt.MaxFunEvals = 2e4;
        P.options(2).nlpOpt.MaxIter = 1e5;
        
    case 'rungeKutta'
        P.options(1).method = 'rungeKutta';
        P.options(1).defaultAccuracy = 'medium';
        
        P.options(2).method = 'rungeKutta';
        P.options(2).defaultAccuracy = 'medium';
        
    case 'hermiteSimpson'
        P.options(1).method = 'hermiteSimpson';
        P.options(1).defaultAccuracy = 'medium';
        
        P.options(2).method = 'hermiteSimpson';
        P.options(2).defaultAccuracy = 'medium';
        
    case 'chebyshev'
        
        P.options(1).method = 'chebyshev';
        P.options(1).defaultAccuracy = 'low';
        
        P.options(2).method = 'chebyshev';
        P.options(2).defaultAccuracy = 'low';
        P.options(2).chebyshev.nColPts = 15;
        
end


%%%% NOTES:
%
% 1) Orthogonal collocation (chebyshev) is not a good method for this problem, beause there is a
% discontinuity in solution of the thrust curve. It still sort of works,
% but will find a sub-optimal answer, or produce ringing.
%
% 2) Why does the 'trapezoid' low resolution version finish so quickly and the medium
% quality one take forever? Hint: Look at the feasibility printout: it is
% cyclical. If you were to plot the solution as a function of iteration,
% you would find that occasionally the discontinuity moves, which causes a
% consistency error in the NLP. Eventually it gets to the "right" answer,
% although it is pretty boring. I suspect that you could get more
% interesting behavior with different constants.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              Solve!                                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soln = optimTraj(P);

t = linspace(soln(end).grid.time(1),soln(end).grid.time(end),250);
x = soln(end).interp.state(t);
u = soln(end).interp.control(t);


%% Plot!! 
% figure(120);
% subplot(2,2,1);
% plot(t,x(1,:))
% xlabel('time (s)')
% ylabel('height [m]')
% grid on
% title('VTOL Take off Trajectory')
% subplot(2,2,2);
% plot(t,x(3,:)*(1-tau))
% xlabel('time (s)')
% ylabel('Energy [J]')
% grid on
% title('VTOL Energy discharge')
% subplot(2,2,3);
% plot(t,x(2,:))
% xlabel('time (s)')
% ylabel('velocity [m/s]')
% grid on
% title('Aircraft velocity distribution')
% subplot(2,2,4);
% plot(t,u)
% xlabel('time (s)')
% ylabel('Engine revolutions [rps]')
% grid on
% title('Engine control')


%% Post process
%% Aerodynamic model
alpha_V = -90*pi/180; % [º] Because of Vertical TO


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


% Full model w.r.t angle of attack

CD = @(alpha) p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + p_CD_ac(4)*alpha.^3 ...
    + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);

%% Thrust
CT0 =  0.0735880531010883;
CT1 = -0.0311758018412727;
CT2 = -0.249744726429543;

CT = @(V,n) CT2*V.^2./(n.^2*D^2)+CT1*V./(n*D)+CT0;

T = @(V,n) CT(V,n).*N_eng*rho.*n.^2*D^4;

CP0 = 0.0261518307541734;
CP1 = 0.0473735972985378;
CP2 = -0.16267474946046;
CP3 = 0.0247028469343899;

CP = @(V,n) CP3*V.^3./(n.^3*D^3)+CP2*V.^2./(n.^2*D^2)+CP1*V./(n*D)+CP0;

P = @(V,n) CP(V,n).*N_eng*rho.*n.^3*D^5;


T_vec = T(x(2,:),u);
P_vec = P(x(2,:),u);
D_vec = 0.5*rho.*x(2,:).^2*S*CD(alpha_V);


%% Plot!! 


figure(1)
plot(t,x(1,:),'LineWidth',1.5)
Title1 = strcat('VTOL Take off Trajectory');
grid on
title(Title1)
xlabel('t [s]')
ylabel('h [m]')

figure(2)
plot(t,x(2,:),'LineWidth',1.5)
Title1 = strcat('Aircraft velocity distribution');
grid on
title(Title1)
xlabel('t [s]')
ylabel('V [m/s]')

figure(3)
plot(t,x(3,:),'LineWidth',1.5)
Title1 = strcat('VTOL Energy discharge');
grid on
title(Title1)
xlabel('t [s]')
ylabel('E [J]')

figure(4)
plot(t,u,'LineWidth',1.5)
Title1 = strcat('Engine control');
grid on
title(Title1)
xlabel('t [s]')
ylabel('Engine revolutions [rps]')

figure(5)
plot(t,T_vec,'LineWidth',1.5)
Title1 = strcat('VTOL Thrust');
grid on
title(Title1)
xlabel('t [s]')
ylabel('T [N]')

figure(6)
plot(t,P_vec,'LineWidth',1.5)
Title1 = strcat('VTOL Power');
grid on
title(Title1)
xlabel('t [s]')
ylabel('P [W]')

figure(7)
plot(t,D_vec,'LineWidth',1.5)
Title1 = strcat('VTOL Aerodynamic Drag');
grid on
title(Title1)
xlabel('t [s]')
ylabel('D [N]')

figure(8)
plot(t,(x(3,:))/E_batt*100,'LineWidth',1.5)
Title1 = strcat('VTOL Energy perc discharge');
grid on
title(Title1)
xlabel('t [s]')
ylabel('E [%]')



%%%%%%%%%%%%%%%%%%%%%
% Extract mission file


t_vec = t;
x_vec = zeros(1,length(t));
h_vec = x(1,:);
V_vec = x(2,:);
n_vec = u;
gamma_vec = pi/2*ones(1,length(t));
alpha_vec = -pi/2*ones(1,length(t));
epsilon_vec = pi/2*ones(1,length(t));
phi_vec = zeros(1,length(t));
L_vec = NaN(1,length(t));
%T_vec = T_vec;
%D_vec = D_vec;
%P_vec = P_vec;
E_vec =  x(3,:);
E_vec_perc = (E_vec)/E_batt*100;

Mission_2 = [t_vec; x_vec ; h_vec; V_vec; n_vec;gamma_vec; alpha_vec ; epsilon_vec; phi_vec; L_vec; T_vec; D_vec; P_vec; E_vec; E_vec_perc]';
save Mission_2.mat Mission_2
