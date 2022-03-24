%% Code for solving Max Range for a given Energy



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


S_ref = 0.430279179101; % Reference surface [m^2]
rho   = 1.225;

tau = 0.2; 
CD0 = 0.035682195723975;
CD2 = 0.054209627025009;

%% Auxiliary variables

alpha_s = N_eng*CT0*D^4;

beta_s = N_eng*CT1*D^3;

psi_s = CT2*D^2-S_ref*CD0/2;

epsilon_s = 2*CD2*(mTOW*g)^2/(rho^2*S_ref);

%% Function Definition

dndV = @(X) -(4*X(1).^3*psi_s + 2.*X(2).^2.*X(1)*alpha_s + 3*beta_s*X(2).*X(1).^2 )./(2*X(2).*X(1).^2*alpha_s + beta_s*X(1).^3);

CP    = @(X) CP0 + CP1*X(1)./(X(2)*D) + CP2*(X(1)./(X(2)*D)).^2 +CP3*(X(1)./(X(2)*D)).^3;
dCPdV = @(X) CP1./(X(2)*D) + 2*CP2*X(1)./(X(2)*D).^2 +3*CP3*X(1).^2./(X(2)*D).^3;
CT    = @(X) CT0 + CT1*X(1)./(X(2)*D) + CT2*(X(1)./(X(2)*D)).^2;
f1    = @(X) CP(X).*X(2) - X(1).*dCPdV(X).*X(2) - 3*CP(X).*dndV(X).*X(1);
f2    = @(X) N_eng*CT(X)-1/2*X(1).^2*S_ref*(CD0+CD2*(2*mTOW*g./(rho*S_ref*X(1))).^2)./(D^4*X(2).^2);
%f5    = @(X) N_eng*(X(2).^2.*X(1).^2*CT0 + CT1*X(1).^3.*X(2)/D + CT2*X(1).^4/D^2) - S_ref/2*X(1).^4*CD0/D^4 - 2*CD2*(mTOW*g)^2/(rho^2*S_ref*D^4);

%%
% X(1) = V
% X(2) = n

options = optimset('TolX',1e-4,'MaxIter',2000);

fsol = @(X)[f1(X);
    f2(X)];

X0 = [30,70];
[Xsol,err] = fsolve(@(X) fsol(X),X0,options)


