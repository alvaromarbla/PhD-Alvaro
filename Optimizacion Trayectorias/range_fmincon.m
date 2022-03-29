%% Code for solving Max Range for a given Energy
clear all
close all
mbatt = 7; %[kg]
e0    = 720e3; %[J/kg]

E     = e0*mbatt;  %[J] Total energy of the battery packs

N_eng = 2; % Number of engines
D     = 0.7112; % Propeller Diameter [m]


CT0 =  0.089050757500461;
CT1 = -0.026659496766654;
CT2 = -0.162620845549081; 


CP0 = 0.034214580122684;
CP1 = -0.002523708791417;
CP2 = 0.116121898742278;
CP3 = -0.248807360672063;

mTOW = 21.39; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.430279179101; % Reference surface [m^2]
rho   = 1.225;

tau = 0.2; 
CD0 = 0.035682195723975;
CD2 = 0.054209627025009;

eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%X(1) = V
%X(2) = n

%% Power
CP = @(X) CP3*X(1)^3/(X(2)^3*D^3)+CP2*X(1)^2/(X(2)^2*D^2)+CP1*X(1)/(X(2)*D)+CP0;

P = @(X) CP(X)*N_eng*rho*X(2)^3*D^5;


%% Range

% - because fmincon solves for the minimum of a function!!
xf = @(X) -eta_m*E*X(1)/((1-tau)*P(X))*1e-3;
nonlcon = @constrains_TD;


x0 = [24.08 55.9];
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];


Xsol = fmincon(xf,x0,A,b,Aeq,beq,lb,ub,nonlcon)


function [c,ceq] = constrains_TD(X)

N_eng = 2; % Number of engines
D     = 0.7112; % Propeller Diameter [m]


CT0 =  0.089050757500461;
CT1 = -0.026659496766654;
CT2 = -0.162620845549081; 


mTOW = 21.39; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.430279179101; % Reference surface [m^2]
rho   = 1.225;

CD0 = 0.035682195723975;
CD2 = 0.054209627025009;

c  = [];
ceq = N_eng*rho*(CT2*X(1)^2/(X(2)^2*D^2)+CT1*X(1)/(X(2)*D)+CT0)*X(2)^2*D^4 ... 
    -(1/2)*rho*X(1)^2*S_ref*(4*W^2*CD2/(rho^2*X(1)^4*S_ref^2)+CD0);
end






