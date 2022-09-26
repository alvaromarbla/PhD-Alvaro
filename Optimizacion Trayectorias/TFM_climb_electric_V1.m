%% Code for solving Max Range for a given Energy
clear all
%close all

% This code takes as inputs:

    % Aerodynamic wind tunnel model // polar simplified model from previous
    % Propulsive  wind tunnel model 
    % Weights parameters from last version (June 2022)
    
% And uses three algorithms:
    % Graphic method
    % fmincon method
    
    % Analytical method not available since problem is too complex and
    % variables can't be isolated

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

%% Set deltaH (altitude difference)

deltaH = 100; %[m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% We start with graphic method 


%% Power
CP = @(V,n) CP3*V^3/(n^3*D^3)+CP2*V^2/(n^2*D^2)+CP1*V/(n*D)+CP0;

P = @(V,n) CP(V,n)*N_eng*rho*n^3*D^5;

%% Thrust
CT = @(V,n) CT2*V^2/(n^2*D^2)+CT1*V/(n*D)+CT0;

T = @(V,n) CT(V,n)*N_eng*rho*n^2*D^4;

%% Lift 

L = @(V,alpha) 0.5*rho*V^2*S_ref*CL(alpha);

%% Drag 

Dr = @(V,alpha) 0.5*rho*V^2*S_ref*CD(alpha);

%% Objective Function: Energy consumption
% gmm = gamma
E_climb = @(V,n,gmm)  P(V,n)*deltaH/(eta_m*(1-tau)*V*sin(gmm))*1e-3;

%% Generate numerical Matrix
Nmat = 250;
Vlin = linspace(15,40,Nmat);
gmmlin = linspace(2*pi/180,75*pi/180, Nmat); % Can't start at 0 because it gives non-defined functions

% Preallocate matrix for speed
Emat = zeros(Nmat,Nmat);
nmat = zeros(Nmat,Nmat);

%% Set fzero for dynamic constraints

options = optimset('TolX',1e-5,'MaxIter',1000);
% Loop in speed
for  ii = 1: Nmat
    
    Xg0 = [60,7*pi/180]; % Just some random but logic data
    ii
    % Loop in gammas
    for jj = 1:Nmat
        % X1 = n
        % X2 = alpha
        fun_graph = @(Xg) [T(Vlin(ii),Xg(1))*cos(Xg(2)+gmmlin(jj))-L(Vlin(ii),Xg(2))*sin(gmmlin(jj))-Dr(Vlin(ii),Xg(2))*cos(gmmlin(jj));...
        T(Vlin(ii),Xg(1))*sin(Xg(2)+gmmlin(jj))+L(Vlin(ii),Xg(2))*cos(gmmlin(jj))-Dr(Vlin(ii),Xg(2))*sin(gmmlin(jj))-W];
        
        % Find rpms (n) and angle of attack (alpha) that solve eqs
        
        Xg_sol = fsolve(@(Xg)fun_graph(Xg),Xg0,options);
        Emat(ii,jj) = E_climb(Vlin(jj),Xg_sol(1),gmmlin(ii));
        nmat(ii,jj) = Xg_sol(1);
        
        Xg0   = Xg_sol;
        
    end
end

%% Finally, solve using fmincon

%% Power
CP_f = @(X) CP3*X(1)^3/(X(2)^3*D^3)+CP2*X(1)^2/(X(2)^2*D^2)+CP1*X(1)/(X(2)*D)+CP0;

P_f = @(X) CP_f(X)*N_eng*rho*X(2)^3*D^5;

%% Energy 
% X1 = V
% X2 = n
% X3 = gamma
% + because fmincon solves for the minimum of a function!!
E_fmincon = @(X) P_f(X)*deltaH/(eta_m*(1-tau)*X(1)*sin(X(3)))*1e-3;


nonlcon = @constrains_TD;


x0_f = [Vlin(1)*1.2 65 5*pi/180 7*pi/180 CL(7*pi/180) CD(7*pi/180) ];
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];


Xsol_fmincon     = fmincon(E_fmincon,x0_f,A,b,Aeq,beq,lb,ub,nonlcon)
nsol_fmincon     = Xsol_fmincon(2);
Vsol_fmincon     = Xsol_fmincon(1);
gammasol_fmincon = Xsol_fmincon(3)*180/pi;
alphasol_fmincon = Xsol_fmincon(4)*180/pi;
CLsol_fmincon    = Xsol_fmincon(5);
CDsol_fmincon    = Xsol_fmincon(6);

%% This maximum considers the stall limitation
E_min_fmin       = E_fmincon (Xsol_fmincon);



%%%%%%
% Now everything is solved!! We only need to plot the results
%%%%%%

%% Plot contours
N_contour_lines = 10; % Number of contour lines
vect_cc_E = linspace(min(min(Emat)),max(max(Emat)),N_contour_lines);
vect_cc_n = linspace(min(min(nmat)),max(max(nmat)),N_contour_lines);


 figure(1)
[Ec_c,h_Ecmat_c] = contourf(gmmlin*180/pi,Vlin,Emat,vect_cc_E');
       clabel(Ec_c,h_Ecmat_c)
       colormap(flipud(colormap('gray')))
         grid on
         Title1 = strcat(['Energy consumption [kJ] vs. V [m/s] & traj. angle [deg]. Min Energy is ', num2str(E_min_fmin), ' kJ']);
         title(Title1)
         xlabel('Trajectory angle [deg] ')
         ylabel('V [m/s]')
         sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
         sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
         sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
         sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
         caxis([min(min(Emat)) max(max(Emat))]); %%%%%%%%%%%%%%%
         hold on 
         plot(gammasol_fmincon,Vsol_fmincon,'og','LineWidth',2.5)
         legend('n(V,\gamma)','Optimal solution.','Location','northwest')

  figure(2)
[n_c,h_nmat_c] = contourf(gmmlin*180/pi,Vlin,nmat,vect_cc_n');
       clabel(n_c,h_nmat_c)
       colormap(flipud(colormap('gray')))
         grid on
         Title1 = strcat('Engine revolutions [rps] vs. V [m/s] & traj. angle [deg]');
         title(Title1)
         xlabel('Trajectory angle [deg] ')
         ylabel('V [m/s]')
         sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
         sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
         sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
         sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
         caxis([min(min(nmat)) max(max(nmat))]); %%%%%%%%%%%%%%%
         hold on 
         plot(gammasol_fmincon,Vsol_fmincon,'og','LineWidth',2.5)
         legend('n(V,\gamma)','Optimal solution.','Location','northwest')
        
         
         
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
%% Lift 

L = @(V,alpha) 0.5*rho*V^2*S_ref*CL(alpha);

%% Drag 

Dr = @(V,alpha) 0.5*rho*V^2*S_ref*CD(alpha);
%X(1) = V
%X(2) = n
%X(3) = gamma 
%X(4) = alpha
%X(5) = CD
%X(6) = CL


% [Stall speed condition; Max RPS condition]
c  = [-CL_max_w1_CR_ope + X(6);-rps_max+X(2)]; 

% Dynamic equations
ceq = [T(X(1),X(2))*cos(X(4)+X(3))-L(X(1),X(4))*sin(X(3))-Dr(X(1),X(4))*cos(X(3));...
       T(X(1),X(2))*sin(X(4)+X(3))+L(X(1),X(4))*cos(X(3))-Dr(X(1),X(4))*sin(X(3))-W;
       CL(X(4))-X(6);
       CD(X(4))-X(5)];
end
