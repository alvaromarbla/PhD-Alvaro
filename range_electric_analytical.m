%% Code for solving Max Range for a given Energy
clear all

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



%% Function Definition
%Derivative of xf w.r.t V
fsol = @(V) -32*eta_m*E*V^3*N_eng^2*rho^2*CT0^3*S_ref^3*D/((1-tau)*(-8*CP3*V^6*CT0^3*N_eng^3*S_ref^3*rho^3*D^3/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^3+...
    4*CP2*V^4*CT0^2*N_eng^2*S_ref^2*rho^2*D^2/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^2-2*CP1*V^2*CT0*N_eng*S_ref*rho*D/...
    (CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))+CP0)*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^3)...
    +8*eta_m*E*V^4*N_eng^2*rho^2*CT0^3*S_ref^3*D*(-48*CP3*V^5*CT0^3*N_eng^3*S_ref^3*rho^3*D^3/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^3+24*CP3*V^6*CT0^3*N_eng^3*S_ref^3*rho^3*D^3*(2*CT1*N_eng*S_ref*D*V*rho-(1/2)*...
    (-16*CT0*CT2*D^2*N_eng^2*S_ref^2*V^3*rho^2+4*CT1^2*D^2*N_eng^2*S_ref^2*V^3*rho^2+8*CD0*CT0*N_eng*S_ref^3*V^3*rho^2)/sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))/...
    (CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^4+16*CP2*V^3*CT0^2*N_eng^2*S_ref^2*rho^2*D^2/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2 ...
    +CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^2-8*CP2*V^4*CT0^2*N_eng^2*S_ref^2*rho^2*D^2*(2*CT1*N_eng*S_ref*D*V*rho-(1/2)*(-16*CT0*CT2*D^2*N_eng^2*S_ref^2*V^3*rho^2+4*CT1^2*D^2*N_eng^2*S_ref^2*V^3*rho^2+8*CD0*CT0*N_eng*S_ref^3*V^3*rho^2)...
    /sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^3-4*CP1*V*CT0*N_eng*S_ref*rho*D/...
    (CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))+2*CP1*V^2*CT0*N_eng*S_ref*rho*D*(2*CT1*N_eng*S_ref*D*V*rho-(1/2)*(-16*CT0*CT2*D^2*N_eng^2*S_ref^2*V^3*rho^2+4*CT1^2*D^2*N_eng^2*S_ref^2*V^3*rho^2+8*CD0*CT0*N_eng*S_ref^3*V^3*rho^2)...
    /sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^2)/((1-tau)*(-8*CP3*V^6*CT0^3*N_eng^3*S_ref^3*rho^3*D^3/...
    (CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^3+4*CP2*V^4*CT0^2*N_eng^2*S_ref^2*rho^2*D^2/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^2-2*CP1*V^2*CT0*N_eng*S_ref*rho*D/...
    (CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))+CP0)^2*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^3)+24*eta_m*E*V^4*N_eng^2*rho^2*CT0^3*S_ref^3*D*(2*CT1*N_eng*S_ref*D*V*rho-(1/2)*...
    (-16*CT0*CT2*D^2*N_eng^2*S_ref^2*V^3*rho^2+4*CT1^2*D^2*N_eng^2*S_ref^2*V^3*rho^2+8*CD0*CT0*N_eng*S_ref^3*V^3*rho^2)/sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))/((1-tau)*(-8*CP3*V^6*CT0^3*N_eng^3*S_ref^3*rho^3*D^3/(CT1*N_eng*S_ref*D*V^2*rho-...
    sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^3+4*CP2*V^4*CT0^2*N_eng^2*S_ref^2*rho^2*D^2/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^2-2*CP1*V^2*CT0*N_eng*S_ref*rho*D...
    /(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))+CP0)*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+8*CT0*N_eng*S_ref*W^2*CD2))^4) ;

%% Set fzero
options = optimset('TolX',1e-5,'MaxIter',2000);


X0 = 25;
[Vsol,err] = fsolve(fsol,X0,options)

%% T= D constrain -> n = n(V)

ncon = @(V) -(1/2)*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+...
    8*CT0*N_eng*S_ref*W^2*CD2))/(CT0*N_eng*S_ref*rho*V*D^2);

nsol = ncon(Vsol)


%% Power
CP = @(V,n) CP3*V^3/(n^3*D^3)+CP2*V^2/(n^2*D^2)+CP1*V/(n*D)+CP0;

P = @(V,n) CP(V,n)*N_eng*rho*n^3*D^5;


%% Range

xf = @(V,n) eta_m*E*V/((1-tau)*P(V,n));
% Evaluate solution
xfsol = xf(Vsol,nsol)





