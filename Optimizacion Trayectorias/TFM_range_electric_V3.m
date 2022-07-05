%% Code for solving Max Range for a given Energy
clear all
close all

% This code takes as imputs:

    % Aerodynamic wind tunnel model // polar simplified model from previous
    % Propulsive  wind tunnel model 
    % Weights parameters from last version (June 2022)
    
% And uses three algorithms:
    % Graphic method
    % Analytical method (for polar only)
    % fmincon method

%% Aerodynamic model
p_CL_ac = [0.582;3.345731;-0.642635;-2.187085;0.713766;0.377985;-0.2946314;-0.004890;0.043976];


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


% Full model w.r.t angle of attack
CL = @(alpha) p_CL_ac(9)*alpha.^8+p_CL_ac(8)*alpha.^7 + p_CL_ac(7)*alpha.^6 + p_CL_ac(6)*alpha.^5 + p_CL_ac(5)*alpha.^4 + p_CL_ac(4)*alpha.^3 ...
    + p_CL_ac(3)*alpha.^2 + p_CL_ac(2)*alpha + p_CL_ac(1);

CD = @(alpha) p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + p_CD_ac(4)*alpha.^3 ...
    + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);
% Parabolic corrected polar (only take CL = CL0 + CLalpha*alpha and CD = CD0+K1*CL+K2*CL^2)

CL0      = p_CL_ac(1);
CL_alpha = p_CL_ac(2);

CD2      = p_CD_ac(3)/CL_alpha^2; % CLO
CD1      = p_CD_ac(2)/CL_alpha - 2*CD2*CL0; % K1
CD0      = p_CD_ac(1)-CD1*CL0-CD2*CL0^2; % K2

alphav = linspace(-pi/2,pi/2,180);


CLv = CL(alphav);
CDv = CD(alphav);
% figure(1)
% plot(alphav*180/pi,CLv)
% figure(2)
% plot(alphav*180/pi,CDv)
% 
% figure(3)
% plot(CLv,CDv)

CL_max_w1_CR  = max(CLv); 
Safety_F      = 1.2;
CL_max_w1_CR_ope = CL_max_w1_CR/Safety_F^2; 


%% Propulsive Model

CT0 =  0.089050757500461;
CT1 = -0.026659496766654;
CT2 = -0.162620845549081; 


CP0 = 0.034214580122684;
CP1 = -0.002523708791417;
CP2 = 0.116121898742278;
CP3 = -0.248807360672063;


%% Loads and weight
mTOW = 21.39; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.430279179101; % Reference surface [m^2]
rho   = 1.225;


%% Electric
tau = 0.2; 
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency
mbatt = 7; %[kg]
e0    = 720e3; %[J/kg]
E     = e0*mbatt;  %[J] Total energy of the battery packs
N_eng = 2; % Number of engines
D     = 0.7112; % Propeller Diameter [m]

RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% We start with graphic method 


%% Power
CP = @(V,n) CP3*V^3/(n^3*D^3)+CP2*V^2/(n^2*D^2)+CP1*V/(n*D)+CP0;

P = @(V,n) CP(V,n)*N_eng*rho*n^3*D^5;

%% Thrust
CT = @(V,n) CT2*V^2/(n^2*D^2)+CT1*V/(n*D)+CT0;

T = @(V,n) CT(V,n)*N_eng*rho*n^2*D^4;


%% Range

xf = @(V,n)  eta_m*E*V/((1-tau)*P(V,n))*1e-3;

%% Generate numerical Matrix
Nmat = 1000;
Vlin = linspace(20,40,Nmat);
nlin = linspace(40,100, Nmat);

% Preallocate matrix for speed
xfmat = zeros(Nmat,Nmat);
% Loop in speed
for  ii = 1: Nmat
    % Loop in revolutions
    for jj = 1:Nmat
        
        xfmat(ii,jj) = xf(Vlin(jj),nlin(ii));
        
    end
end

%% Calculate restriction T = D 


% Using simplified aero model

ncon = @(V) -(1/2)*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+...
    CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))/(CT0*N_eng*S_ref*rho*V*D^2);
% Preallocate matrix for speed
nconlin = zeros(Nmat,1);
for kk = 1:Nmat
    nconlin (kk) = ncon(Vlin(kk));

end


% Using full aero model

CLcon = @(V) W./(1/2*rho*V.^2*S_ref); 
X_sol_full = zeros(Nmat,2);
    %X1 = alpha
    %X2 = n
    
    alpha0 = 6*pi/180; % First iteration
    X0_full = [alpha0,nconlin(1)];
    options = optimset('TolX',1e-5,'MaxIter',1000);
   
    for kk = 1:Nmat
        % Define function [ L = W, T = D]
        fun_full = @(X_full)[-CLcon(Vlin(kk))+CL(X_full(1));T(Vlin(kk),X_full(2))-0.5*rho*Vlin(kk)^2*S_ref*CD(X_full(1))];


        
        X_sol_full(kk,:) = fsolve(@(X_full)fun_full(X_full),X0_full,options);
        % Update initial iterant
        X0_full = X_sol_full(kk,:)  ;
    end
    % Extract revolutions for full model
nconlin_full = X_sol_full(:,2);
    

%% Calculate stall speed

Vstall = sqrt(2*W/(rho*S_ref*CL_max_w1_CR_ope));

%% Eliminate zeros from indeterminations (points where CP = 0)
xfmat  = xfmat'; % To make revolutions stand on the x-axis and V on the y-axis


threshold = 2e2; % Assume Range < 200 km 
xfmat(xfmat > threshold) = 0;
xfmat(xfmat < 0) = 0;


 
%% Now, solve in an Analytical manner


%% Function Definition
%Derivative of xf w.r.t V
fsol = @(V) -32*eta_m*E*V^3*N_eng^2*rho^2*CT0^3*S_ref^3*D/((1-tau)*(-8*CP3*V^6*CT0^3*N_eng^3*S_ref^3*rho^3*D^3/(CT1*N_eng*S_ref*D*V^2*rho...
    -sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^3 ...
    +4*CP2*V^4*CT0^2*N_eng^2*S_ref^2*rho^2*D^2/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2 ...
    +2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^2-2*CP1*V^2*CT0*N_eng*S_ref*rho*D/(CT1*N_eng*S_ref*D*V^2*rho... 
    -sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho ...
    +8*CD2*CT0*N_eng*S_ref*W^2))+CP0)*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2 ... 
    +2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^3)+8*eta_m*E*V^4*N_eng^2*rho^2*CT0^3*S_ref^3*D*(-48*CP3*V^5 ...
    *CT0^3*N_eng^3*S_ref^3*rho^3*D^3/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*...
    S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^3+24*CP3*V^6*CT0^3*N_eng^3*S_ref^3*rho^3*D^3*(2*CT1*N_eng*S_ref*D*V*rho-(1/2)*... 
    (-16*CT0*CT2*D^2*N_eng^2*S_ref^2*V^3*rho^2+4*CT1^2*D^2*N_eng^2*S_ref^2*V^3*rho^2+8*CD0*CT0*N_eng*S_ref^3*V^3*rho^2+8*CD1*CT0*N_eng*S_ref^2*V*W*rho)/sqrt(-4*CT0*... 
    CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))...
    /(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*...
    S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^4+16*CP2*V^3*CT0^2*N_eng^2*S_ref^2*rho^2*D^2/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*...
    rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^2-8*CP2*V^4*CT0^2*...
    N_eng^2*S_ref^2*rho^2*D^2*(2*CT1*N_eng*S_ref*D*V*rho-(1/2)*(-16*CT0*CT2*D^2*N_eng^2*S_ref^2*V^3*rho^2+4*CT1^2*D^2*N_eng^2*S_ref^2*V^3*rho^2+8*CD0*CT0*N_eng*...
    S_ref^3*V^3*rho^2+8*CD1*CT0*N_eng*S_ref^2*V*W*rho)/sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4 ... 
    *rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*...
    N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^3-4*CP1*V*CT0*N_eng*S_ref*rho*D/...
    (CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng...
    *S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))+2*CP1*V^2*CT0*N_eng*S_ref*rho*D*(2*CT1*N_eng*S_ref*D*V*rho-(1/2)*(-16*CT0*CT2*D^2*N_eng^2*S_ref^2*V^3*rho^2+...
    4*CT1^2*D^2*N_eng^2*S_ref^2*V^3*rho^2+8*CD0*CT0*N_eng*S_ref^3*V^3*rho^2+8*CD1*CT0*N_eng*S_ref^2*V*W*rho)/sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2 ...
    *D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))/(CT1*N_eng*S_ref*D*V^2*rho-...
    sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+...
    8*CD2*CT0*N_eng*S_ref*W^2))^2)/((1-tau)*(-8*CP3*V^6*CT0^3*N_eng^3*S_ref^3*rho^3*D^3/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4 ...
    *rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^3+...
    4*CP2*V^4*CT0^2*N_eng^2*S_ref^2*rho^2*D^2/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2 ...
    +2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^2-2*CP1*V^2*CT0*N_eng*S_ref*rho*D/(CT1*N_eng*S_ref*D*V^2*rho-...
    sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+...
    8*CD2*CT0*N_eng*S_ref*W^2))+CP0)^2*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0...
    *N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^3)+24*eta_m*E*V^4*N_eng^2*rho^2*CT0^3*S_ref^3*D*(2*CT1*N_eng*S_ref*D*V*rho-...
    (1/2)*(-16*CT0*CT2*D^2*N_eng^2*S_ref^2*V^3*rho^2+4*CT1^2*D^2*N_eng^2*S_ref^2*V^3*rho^2+8*CD0*CT0*N_eng*S_ref^3*V^3*rho^2+8*CD1*CT0*N_eng*S_ref^2*V*W*rho)/...
    sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+...
    8*CD2*CT0*N_eng*S_ref*W^2))/((1-tau)*(-8*CP3*V^6*CT0^3*N_eng^3*S_ref^3*rho^3*D^3/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2 ...
    +CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^3+4*CP2*V^4*CT0^2*N_eng^2 ...
    *S_ref^2*rho^2*D^2/(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2 ...
    +4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^2-2*CP1*V^2*CT0*N_eng*S_ref*rho*D/(CT1*N_eng*S_ref*D*V^2*rho-...
    sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+...
    8*CD2*CT0*N_eng*S_ref*W^2))+CP0)*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*...
    N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))^4);




%% Set fsolve
options = optimset('TolX',1e-5,'MaxIter',2000);


X0 = 25;
[Vsol,err] = fsolve(fsol,X0,options);

ncon_A = @(V) -(1/2)*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+...
    2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+4*CD1*CT0*N_eng*S_ref^2*V^2*W*rho+8*CD2*CT0*N_eng*S_ref*W^2))/(CT0*N_eng*S_ref*rho*V*D^2);

% Calculate the revolutions from the T = D restriction
nsol_A = ncon_A(Vsol);

% Solve this value
xf_sol_A = xf(Vsol,nsol_A);


%% Finally, solve using fmincon

%% Power
CP_f = @(X) CP3*X(1)^3/(X(2)^3*D^3)+CP2*X(1)^2/(X(2)^2*D^2)+CP1*X(1)/(X(2)*D)+CP0;

P_f = @(X) CP_f(X)*N_eng*rho*X(2)^3*D^5;

%% Range

% - because fmincon solves for the minimum of a function!!
xf_fmincon = @(X) -eta_m*E*X(1)/((1-tau)*P_f(X))*1e-3;
nonlcon = @constrains_TD;


x0_f = [Vlin(1) nconlin(1) X_sol_full(1,1) CL(X_sol_full(1,1)) CD(X_sol_full(1,1)) ];
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];


Xsol_fmincon = fmincon(xf_fmincon,x0_f,A,b,Aeq,beq,lb,ub,nonlcon)
nsol_fmincon = Xsol_fmincon(2);
Vsol_fmincon = Xsol_fmincon(1);

%% This maximum considers the stall limitation
xf_max_fmin       = -xf_fmincon (Xsol_fmincon);



%%%%%%
% Now everything is solved!! We only need to plot the results
%%%%%%

%% Plot contours
N_contour_lines = 10; % Number of contour lines
vect_cc_xf = linspace(min(min(xfmat)),max(max(xfmat)),N_contour_lines);


 figure(1)
[xf_c,h_xfmat_c] = contourf(nlin,Vlin,xfmat,vect_cc_xf');
       clabel(xf_c,h_xfmat_c)
       colormap(flipud(colormap('gray')))
         grid on
         Title1 = strcat(['Range [km] vs. V [m/s] & Engine rev. [rps]. Max Range is ', num2str((xf_max_fmin+xf_sol_A)/2), ' km']); % Ponderate both results (A+Fmincon) by 1/2
         title(Title1)
         xlabel('Engine revolutions [rps] ')
         ylabel('V [m/s]')
         sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
         sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
         sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
         sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
         caxis([min(min(xfmat)) max(max(xfmat))]); %%%%%%%%%%%%%%%
         hold on 
         plot(nconlin,Vlin,'--r','LineWidth',2)
         plot(nconlin_full,Vlin,'-.m','LineWidth',2)
         yline(Vstall,'-.b','LineWidth',2)
         xline(rps_max,':b','LineWidth',2)
         plot(nsol_A,Vsol,'*g','LineWidth',2)
         plot(nsol_fmincon,Vsol_fmincon,'og','LineWidth',2.5)
         legend('x_f (V,n)','T = D constrain. App Model','T = D constrain. Full Model','Min. Operative Speed','Max RPS','Optimal solution. App Model','Optimal solution. Full Model','Location','northwest')

% figure(2)
% 
% mesh(Vlin,nlin,xfmat)



%% Auxiliary function for fmincon
function [c,ceq] = constrains_TD(X)
%% Aerodynamic model
p_CL_ac = [0.582;3.345731;-0.642635;-2.187085;0.713766;0.377985;-0.2946314;-0.004890;0.043976];
p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];

CL = @(alpha) p_CL_ac(9)*alpha.^8+p_CL_ac(8)*alpha.^7 + p_CL_ac(7)*alpha.^6 + p_CL_ac(6)*alpha.^5 + p_CL_ac(5)*alpha.^4 + p_CL_ac(4)*alpha.^3 ...
    + p_CL_ac(3)*alpha.^2 + p_CL_ac(2)*alpha + p_CL_ac(1);

CD = @(alpha) p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + p_CD_ac(4)*alpha.^3 ...
    + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);


% Stall calculations


alphav = linspace(-pi/2,pi/2,180);

CLv = CL(alphav);
CL_max_w1_CR  = max(CLv);

Safety_F      = 1.2;
CL_max_w1_CR_ope = CL_max_w1_CR/Safety_F^2;


%%%%%%
N_eng = 2; % Number of engines
D     = 0.7112; % Propeller Diameter [m]
RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter


%% Weight
mTOW = 21.39; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.430279179101; % Reference surface [m^2]
rho   = 1.225;

%% Thrust
CT0 =  0.089050757500461;
CT1 = -0.026659496766654;
CT2 = -0.162620845549081; 

CT = @(V,n) CT2*V^2/(n^2*D^2)+CT1*V/(n*D)+CT0;

T = @(V,n) CT(V,n)*N_eng*rho*n^2*D^4;

%X(1) = V
%X(2) = n
%X(3) = alpha
%X(4) = CD
%X(5) = CL


% [Stall speed condition; Max RPS condition]s
c  = [+sqrt(2*W/(rho*S_ref*CL_max_w1_CR_ope))-X(1);-rps_max+X(2)]; 

% T = D condition
ceq = [ T(X(1),X(2))- 0.5*rho*X(1)^2*S_ref*X(4); 
    sqrt(2*W/(rho*S_ref*X(5)))-X(1);
    CL(X(3))-X(5);
    CD(X(3))-X(4)];
end
