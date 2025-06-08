%% Code for solving Max Range for a given Energy
clear all
close all

% This code takes as inputs:

    % Aerodynamic wind tunnel model // polar simplified model from previous
    % Propulsive  wind tunnel model (Marta Nuñez, 2022). 2 dof model
    % collapsed to phi = 0.
    
    % Weights parameters from last version (June 2022)

  % fmincon method




%% Loads and weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.5008; % Reference surface [m^2]
rho   = 1.2133;


%% Electric
tau = 0.0; 
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency
mbatt = 3; %[kg]
e0    = 720e3; %[J/kg]
E     = e0*mbatt*0.6  %[J] Total energy of the battery packs
%E = 1.2111e+06 - e0*mbatt*0.2; % Total energy after mission to perform extended cruise
N_eng = 2; % Number of engines
D     = 0.8128*0.85; % Propeller Diameter [m]

RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;

%% Finally, solve using fmincon

%% Power
%% Propulsive Model

a1_P = 0.0261518307541734;
b1_P = 0.0473735972985378;
c1_P = -0.16267474946046;
d1_P = 0.0247028469343899;
e1_P = 0.0306053713439883;

b2_P = -0.0762350484603968;
c2_P = 0.148580471912353;
d2_P = -0.0726017200715775;

b3_P = 0.0897273366920878;
c3_P = 0.0122602815262456;

b4_P = -0.0486029866039398;

J   = @(V,n) V./(n*D);

phi = @(alpha, eps_eng) alpha + eps_eng;

CP_f = @(V,n,alpha,eps_eng)  a1_P + b1_P*J(V,n) + c1_P*J(V,n).^2 + d1_P*J(V,n).^3 + e1_P*J(V,n).^4 + ...
    phi(alpha,eps_eng).*(b2_P*J(V,n) + c2_P*J(V,n).^2 + d2_P*J(V,n).^3) +  phi(alpha,eps_eng).^2.*(b3_P*J(V,n) + c3_P*J(V,n).^2)...
    +  phi(alpha,eps_eng).^3.*b4_P*J(V,n);

P_f = @(X) CP_f(X(1),X(2),X(3),X(6))*N_eng*rho*X(2)^3*D^5;

%% Range

% - because fmincon solves for the minimum of a function!!
xf_fmincon = @(X) -eta_m*E*X(1)/(P_f(X))*1e-3;
nonlcon = @constrains_TD;

%[V n alpha 0.089094147 0.73578832 eps]
x0_f = [24.225109 25.214554 0.046444 0.089094147 0.73578832 1.0593];
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];

%X(1) = V
%X(2) = n
%X(3) = alpha
%X(4) = CD
%X(5) = CL
%X(6) = eps_eng

Xsol_fmincon = fmincon(xf_fmincon,x0_f,A,b,Aeq,beq,lb,ub,nonlcon)
nsol_fmincon = Xsol_fmincon(2);
Vsol_fmincon = Xsol_fmincon(1);
alphasol_fmincon = Xsol_fmincon(3);
CDsol_fmincon = Xsol_fmincon(4);
CLsol_fmincon = Xsol_fmincon(5);
epsilonsol_fmincon = Xsol_fmincon(6);

Lsol = 0.5*rho*Vsol_fmincon^2*S_ref*CLsol_fmincon;
Dsol = 0.5*rho*Vsol_fmincon^2*S_ref*CDsol_fmincon;


%% This maximum considers the stall limitation
xf_max_fmin       = -xf_fmincon (Xsol_fmincon);



%%%%%%%%%%%%%%%%%%%%% Prepare for plots

%% Generate numerical Matrix
Nmat = 100;
Vlin = linspace(17,35,Nmat);
nlin = linspace(20,40,Nmat);
alphalin = linspace(-5,45,Nmat)*pi/180;
epsilonlin = linspace(-5,85, Nmat)*pi/180;

P_f_plot = @(V,n,alpha,eps_eng) CP_f(V,n,alpha,eps_eng)*N_eng*rho*n^3*D^5;
xf_fmincon_plot = @(V,n,alpha,eps_eng) eta_m*E*V/(P_f_plot(V,n,alpha,eps_eng))*1e-3;




% Plot n vs V (1)

% Preallocate matrix for speed
xfmat1 = zeros(Nmat,Nmat);
% Loop in speed
for  ii = 1: Nmat
    % Loop in revolutions
    for jj = 1:Nmat
        
        xfmat1(ii,jj) = xf_fmincon_plot(Vlin(jj),nlin(ii),alphasol_fmincon,epsilonsol_fmincon);
        
    end
end

xf_max_fmin
% %% Eliminate zeros from indeterminations (points where CP = 0)
% xfmat1  = xfmat1'; % To make revolutions stand on the x-axis and V on the y-axis
% 
% 
% threshold = 4.0e2; % Assume Range < 300 km 
% xfmat1(xfmat1 > threshold) = 0;
% xfmat1(xfmat1 < 0) = 0;
% 
% %% Calculate restrictions
% 
% [Vplot1,nplot1] = calculate_restrictions1 (Nmat, Vlin,nlin,Xsol_fmincon);

% Hor forces 1
% Ver forces 2

% Set ratios and figure properties
% Create figure
% hFIG1 = figure;
% fname = "range_electric_V5_1";
% 
% vect_cc_xf = [30,40,50,60,80,80,xf_max_fmin,100,120,140,200]; % Number of contour lines
% [xf_c,h_xfmat1_c] = contourf(nlin,Vlin,xfmat1,vect_cc_xf');
%        clabel(xf_c,h_xfmat1_c)
%        colormap(flipud(colormap('gray')))
%          grid on
%          Title1 = strcat(['Range [km] vs. V [m/s] \& Engine rev. [rps]. ']); % Ponderate both results (A+Fmincon) by 1/2
%          title(Title1)
%          xlabel('Engine revolutions [rps] ')
%          ylabel('V [m/s]')
%          sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
%          sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
%          sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
%          sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
%          %caxis([min(min(xfmat1)) max(max(xfmat1))]); %%%%%%%%%%%%%%%
%          axis([min(nlin) max(nlin) min(Vlin) max(Vlin)])
%          hold on 
%          plot(nplot1(1,:),Vplot1(1,:),'--r','LineWidth',2)
%          plot(nplot1(2,:),Vplot1(2,:),'-.m','LineWidth',2)
%          %plot(nplot(3,:),Vplot(3,:),'-b','LineWidth',2.5)
%          plot(nsol_fmincon,Vsol_fmincon,'*g','LineWidth',2)
%          legend('$x_f$ (V, n)','Hor. Forces cons.','Ver. Forces cons.','Optimal solution','Location','northeast',"Box","on")
% 
% 
%  % Create title
% picturewidth = 30;
% hw_ratio = 0.5;
% 
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex')
% set(findall(hFIG1,'-property','FontSize'),'FontSize',16)
% set(findall(hFIG1,'-property','Box'),'Box','on') % optional
% set(findall(hFIG1,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hFIG1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(hFIG1,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(hFIG1,'Position');
% set(hFIG1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% print(hFIG1,fname,'-dpdf','-vector','-fillpage')        
% 
% % Plot n vs eps (2)legend('x_f (V,n)','Hor. Forces constrain','Ver. Forces constrain','Optimal solution','Location','northwest')
% 
% % Preallocate matrix for speed
% xfmat2 = zeros(Nmat,Nmat);
% % Loop in speed
% for  ii = 1: Nmat
%     % Loop in revolutions
%     for jj = 1:Nmat
% 
%         xfmat2(ii,jj) = xf_fmincon_plot(Vsol_fmincon,nlin(ii),alphasol_fmincon,epsilonlin(jj));
% 
%     end
% end
% 
% %% Eliminate zeros from indeterminations (points where CP = 0)
% xfmat2  = xfmat2'; % To make revolutions stand on the x-axis and V on the y-axis
% epsilonlindeg = epsilonlin*180/pi;
% 
% threshold = 4.0e3; % Assume Range < 300 km 
% xfmat2(xfmat2 > threshold) = 0;
% xfmat2(xfmat2 < 0) = 0;
% 
% %% Calculate restrictions
% 
% 
% [epsplot2,nplot2] = calculate_restrictions2 (Nmat, Vlin,nlin,Xsol_fmincon);
% 
% % Hor forces 1
% % Ver forces 2
% 
% % Set ratios and figure properties
% % Create figure
% hFIG2 = figure;
% fname = "range_electric_V5_2";
% 
% vect_cc_xf = [30,40,50,60,80,80,xf_max_fmin,100,120,140,200]; % Number of contour lines
% [xf_c,h_xfmat2_c] = contourf(nlin,epsilonlindeg,xfmat2,vect_cc_xf');
%        clabel(xf_c,h_xfmat2_c)
%        colormap(flipud(colormap('gray')))
%          grid on
%          Title1 = strcat(['Range [km] vs. $\varepsilon$ [deg] \& Engine rev. [rps]. ']); % Ponderate both results (A+Fmincon) by 1/2
%          title(Title1)
%          xlabel('Engine revolutions [rps] ')
%          ylabel('Engine tilt angle [deg]')
%          sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
%          sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
%          sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
%          sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
%          %caxis([min(min(xfmat2)) max(max(xfmat2))]); %%%%%%%%%%%%%%%
%          hold on 
%          axis([min(nlin) max(nlin) min(epsilonlin*180/pi) max(epsilonlin*180/pi)])
%          hold on 
%          plot(nplot2(1,:),epsplot2(1,:)*180/pi,'--r','LineWidth',2)
%          plot(nplot2(2,:),epsplot2(2,:)*180/pi,'-.m','LineWidth',2)
%          %plot(nplot(3,:),Vplot(3,:),'-b','LineWidth',2.5)
%          plot(nsol_fmincon,epsilonsol_fmincon*180/pi,'*g','LineWidth',2)
%          legend('$x_f$ ($\varepsilon$,n)','Hor. Forces cons.','Ver. Forces cons.','Optimal solution','Location','northeast',"Box","on")
% 
% % Create title
% picturewidth = 30;
% hw_ratio = 0.5;
% 
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex')
% set(findall(hFIG2,'-property','FontSize'),'FontSize',16)
% set(findall(hFIG2,'-property','Box'),'Box','on') % optional
% set(findall(hFIG2,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hFIG2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(hFIG2,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(hFIG2,'Position');
% set(hFIG2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% print(hFIG2,fname,'-dpdf','-vector','-fillpage')
% % Plot V vs alpha (3)
% 
% % Preallocate matrix for speed
% xfmat3 = zeros(Nmat,Nmat);
% % Loop in speed
% for  ii = 1: Nmat
%     % Loop in revolutions
%     for jj = 1:Nmat
% 
%         xfmat3(ii,jj) = xf_fmincon_plot(Vlin(jj),nsol_fmincon,alphalin(ii),epsilonsol_fmincon);
% 
%     end
% end
% 
% %% Eliminate zeros from indeterminations (points where CP = 0)
% xfmat3  = xfmat3'; % To make revolutions stand on the x-axis and V on the y-axis
% alphalindeg = alphalin*180/pi;
% 
% threshold = 4.0e3; % Assume Range < 300 km 
% xfmat3(xfmat3 > threshold) = 0;
% xfmat3(xfmat3 < 0) = 0;
% 
% %% Calculate restrictions
% 
% [alphaplot3,Vplot3] = calculate_restrictions3 (Nmat, alphalin,Vlin,Xsol_fmincon);
% 
% 
% 
% % Set ratios and figure properties
% % Create figure
% hFIG3 = figure;
% fname = "range_electric_V5_3";
% 
% vect_cc_xf = [30,40,50,60,80,80,xf_max_fmin,100,120,140,200]; % Number of contour lines
% [xf_c,h_xfmat3_c] = contourf(alphalindeg,Vlin,xfmat3,vect_cc_xf');
%        clabel(xf_c,h_xfmat3_c)
%        colormap(flipud(colormap('gray')))
%          grid on
%          Title1 = strcat(['Range [km] vs. V [m/s] \&  $\alpha$ [deg]. ']); % Ponderate both results (A+Fmincon) by 1/2
%          title(Title1)
%          xlabel('Angle of attack [deg] ')
%          ylabel('V [m/s]')
%          sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
%          sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
%          sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
%          sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
%          %caxis([min(min(xfmat3)) max(max(xfmat3))]); %%%%%%%%%%%%%%%
%          hold on 
%          axis([min(alphalin*180/pi) max(alphalin*180/pi) min(Vlin) max(Vlin)])
%          hold on 
%          plot(alphaplot3(1,:)*180/pi,Vplot3(1,:),'--r','LineWidth',2)
%          plot(alphaplot3(2,:)*180/pi,Vplot3(2,:),'-.m','LineWidth',2)
%          %plot(nplot(3,:),Vplot(3,:),'-b','LineWidth',2.5)
%          plot(alphasol_fmincon*180/pi,Vsol_fmincon,'*g','LineWidth',2)
%          legend('$x_f$ (V,$\alpha$)','Hor. Forces cons.','Ver. Forces cons.','Optimal solution','Location','northeast',"Box","on")
% 
% % Create title
% picturewidth = 30;
% hw_ratio = 0.5;
% 
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex')
% set(findall(hFIG3,'-property','FontSize'),'FontSize',16)
% set(findall(hFIG3,'-property','Box'),'Box','on') % optional
% set(findall(hFIG3,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hFIG3,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(hFIG3,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(hFIG3,'Position');
% set(hFIG3,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% print(hFIG3,fname,'-dpdf','-vector','-fillpage')
%  % Plot eps vs alpha (4)
% 
% % Preallocate matrix for speed
% xfmat4 = zeros(Nmat,Nmat);
% % Loop in speed
% for  ii = 1: Nmat
%     % Loop in revolutions
%     for jj = 1:Nmat
% 
%         xfmat4(ii,jj) = xf_fmincon_plot(Vsol_fmincon,nsol_fmincon,alphalin(ii),epsilonlin(jj));
% 
%     end
% end
% 
% 
% %% Eliminate zeros from indeterminations (points where CP = 0)
% xfmat4  = xfmat4'; % To make revolutions stand on the x-axis and V on the y-axis
% alphalindeg = alphalin*180/pi;
% 
% threshold = 4.0e3; % Assume Range < 300 km 
% xfmat4(xfmat4 > threshold) = 0;
% xfmat4(xfmat4 < 0) = 0;
% 
% %% Calculate restrictions
% 
% [alphaplot4,epsplot4] = calculate_restrictions4 (Nmat, alphalin,epsilonlin,Xsol_fmincon);
% 
% % Set ratios and figure properties
% % Create figure
% hFIG4 = figure;
% fname = "range_electric_V5_4";
% 
% 
% vect_cc_xf = [30,40,50,60,80,80,xf_max_fmin,100,120,140,200]; % Number of contour lines
% [xf_c,h_xfmat4_c] = contourf(alphalindeg,epsilonlindeg,xfmat4,vect_cc_xf');
%        clabel(xf_c,h_xfmat4_c)
%        colormap(flipud(colormap('gray')))
%          grid on
%          Title1 = strcat(['Range [km] vs. $\varepsilon$ [deg] \&  $\alpha$ [deg].']); % Ponderate both results (A+Fmincon) by 1/2
%          title(Title1)
%          xlabel('Angle of attack [deg] ')
%          ylabel('Engine tilt angle [deg]')
%          sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
%          sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
%          sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
%          sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
%          %caxis([min(min(xfmat4)) max(max(xfmat4))]); %%%%%%%%%%%%%%%
%          hold on
%          axis([min(alphalin*180/pi) max(alphalin*180/pi) min(epsilonlin)*180/pi max(epsilonlin)*180/pi])
%          hold on 
%          plot(alphaplot4(1,:)*180/pi,epsplot4(1,:)*180/pi,'--r','LineWidth',2)
%          plot(alphaplot4(2,:)*180/pi,epsplot4(2,:)*180/pi,'-.m','LineWidth',2)
%          %plot(nplot(3,:),Vplot(3,:),'-b','LineWidth',2.5)
%          plot(alphasol_fmincon*180/pi,epsilonsol_fmincon*180/pi,'*g','LineWidth',2)
%          legend('$x_f$ ($\alpha$,$\varepsilon$)','Hor. Forces cons.','Ver. Forces cons.','Optimal solution','Location','northeast',"Box","on")
% 
% %Create title
% %sgtitle([". Max Range is ", num2str(xf_max_fmin), ' km']);
% % Create title
% picturewidth = 30;
% hw_ratio = 0.5;
% 
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex')
% set(findall(hFIG4,'-property','FontSize'),'FontSize',14)
% set(findall(hFIG4,'-property','Box'),'Box','on') % optional
% set(findall(hFIG4,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hFIG4,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(hFIG4,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(hFIG4,'Position');
% set(hFIG4,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% print(hFIG4,fname,'-dpdf','-vector','-fillpage')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
D     = 0.8128; % Propeller Diameter [m]
RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;

%% Weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.5008; % Reference surface [m^2]
rho   = 1.2133;

%% Thrust


%% Propulsive Model

J   = @(V,n) V./(n*D);
phi = @(alpha, eps_eng) alpha + eps_eng;

 a1_T = 0.0735880531010883;
    b1_T = -0.0311758018412727;
    c1_T = -0.249744726429543;
    d1_T = 0.143084420694372;
    e1_T = 0.0261032283758581;

    b2_T = -0.0982459868751664;
    c2_T = 0.20127470719351;
    d2_T = -0.173738749783189;

    b3_T = 0.156239779501715;
    c3_T = 0.0368592048084175;

    b4_T = -0.0478709034281346;

    CT = @(V,n,alpha,eps_eng)  a1_T + b1_T*J(V,n) + c1_T*J(V,n).^2 + d1_T*J(V,n).^3 + e1_T*J(V,n).^4 + ...
        phi(alpha,eps_eng).*(b2_T*J(V,n) + c2_T*J(V,n).^2 + d2_T*J(V,n).^3) +  phi(alpha,eps_eng).^2.*(b3_T*J(V,n) + c3_T*J(V,n).^2)...
        +  phi(alpha,eps_eng).^3.*b4_T*J(V,n);

    T = @(V,n,alpha,eps_eng) CT(V,n,alpha,eps_eng).*N_eng*rho*n.^2*D^4;

%X(1) = V
%X(2) = n
%X(3) = alpha
%X(4) = CD
%X(5) = CL
%X(6) = eps_eng


% [Stall speed condition; Max RPS condition]s
c  = [+sqrt(2*W/(rho*S_ref*CL_max_w1_CR_ope))-X(1);-rps_max+X(2);T(X(1),X(2),X(3),X(6))-Tmax;X(3)-pi/2;-X(4);-X(6);X(6)-pi/2]; 

% T = D condition
ceq = [ T(X(1),X(2),X(3),X(6)).*cos(phi(X(3),X(6)))- 0.5*rho*X(1)^2*S_ref*X(4); 
    T(X(1),X(2),X(3),X(6)).*sin(phi(X(3),X(6)))+0.5*rho*X(1).^2*S_ref*X(5)-W;
    CL(X(3))-X(5);
    CD(X(3))-X(4)];
end

function [Vplot1,nplot1] = calculate_restrictions1 (Nmat, Vlin,nlin,Xsol_fmincon)
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
D     = 0.8128; % Propeller Diameter [m]
RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;

%% Weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.5008; % Reference surface [m^2]
rho   = 1.2133;

%% Thrust


%% Propulsive Model

J   = @(V,n) V./(n*D);
phi = @(alpha, eps_eng) alpha + eps_eng;

 a1_T = 0.0735880531010883;
    b1_T = -0.0311758018412727;
    c1_T = -0.249744726429543;
    d1_T = 0.143084420694372;
    e1_T = 0.0261032283758581;

    b2_T = -0.0982459868751664;
    c2_T = 0.20127470719351;
    d2_T = -0.173738749783189;

    b3_T = 0.156239779501715;
    c3_T = 0.0368592048084175;

    b4_T = -0.0478709034281346;

    CT = @(V,n,alpha,eps_eng)  a1_T + b1_T*J(V,n) + c1_T*J(V,n).^2 + d1_T*J(V,n).^3 + e1_T*J(V,n).^4 + ...
        phi(alpha,eps_eng).*(b2_T*J(V,n) + c2_T*J(V,n).^2 + d2_T*J(V,n).^3) +  phi(alpha,eps_eng).^2.*(b3_T*J(V,n) + c3_T*J(V,n).^2)...
        +  phi(alpha,eps_eng).^3.*b4_T*J(V,n);

    T = @(V,n,alpha,eps_eng) CT(V,n,alpha,eps_eng).*N_eng*rho*n.^2*D^4;

%X(1) = V
%X(2) = n
%X(3) = alpha
%X(4) = CD
%X(5) = CL
%X(6) = eps_eng

Vsol_fmincon = Xsol_fmincon(1);
nsol_fmincon = Xsol_fmincon(2);
alphasol_fmincon = Xsol_fmincon(3);
CLsol_fmincon = Xsol_fmincon(5);
CDsol_fmincon = Xsol_fmincon(4);
epsilonsol_fmincon = Xsol_fmincon(6);
% Constraint 1 : Hor forces

Vcons1 = zeros(1,Nmat);
ncons1 = nlin;

%func1 = @(V) [ T(V,ncons1(jj),alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*V^2*S_ref*CD(alphasol_fmincon); 
   % T(X(1),X(2),alphasol_fmincon,epsilonsol_fmincon).*sin(phi(alphasol_fmincon,epsilonsol_fmincon))+0.5*rho*X(1).^2*S_ref*CL(alphasol_fmincon)-W];

for jj = Nmat:-1:1
    if jj == Nmat
        x0 = Vsol_fmincon;
    else
        x0 = Vcons1(jj+1);
    end
Vcons1(Nmat-jj+1)= fsolve(@(V) T(V,ncons1(Nmat-jj+1),alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*V.^2*S_ref*CD(alphasol_fmincon), x0);
end
% disp("Valor falso?")
% T(24.20,14.54,alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*24.20.^2*S_ref*CD(alphasol_fmincon)
% disp("Valor verdadero?")
% T(24.2251,25.2146,alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*Vsol_fmincon.^2*S_ref*CD(alphasol_fmincon)
% Constraint 2: Ver Forces
Vcons2 = zeros(1,Nmat);
ncons2 = nlin;


for jj = 1:Nmat
    if jj == 1
        x0 = Vsol_fmincon;
    else
        x0 = Vcons2(jj-1);
    end
Vcons2(jj)= fsolve(@(V) T(V,ncons2(jj),alphasol_fmincon,epsilonsol_fmincon).*sin(phi(alphasol_fmincon,epsilonsol_fmincon))+0.5*rho*V.^2*S_ref*CL(alphasol_fmincon)-W,x0);
end


% Constraint 3: T = Tmax

Vcons3 = ones(1,Nmat)*sqrt(Tmax.*cos(phi(alphasol_fmincon,epsilonsol_fmincon))/(0.5*rho*S_ref*CD(alphasol_fmincon)));
ncons3 = nlin;



Vplot1 = [Vcons1; Vcons2;Vcons3];
nplot1 = [ncons1; ncons2;ncons3];


end

function [epsplot2,nplot2] = calculate_restrictions2 (Nmat, Vlin,nlin,Xsol_fmincon)
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
D     = 0.8128; % Propeller Diameter [m]
RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;

%% Weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.5008; % Reference surface [m^2]
rho   = 1.2133;

%% Thrust


%% Propulsive Model

J   = @(V,n) V./(n*D);
phi = @(alpha, eps_eng) alpha + eps_eng;

 a1_T = 0.0735880531010883;
    b1_T = -0.0311758018412727;
    c1_T = -0.249744726429543;
    d1_T = 0.143084420694372;
    e1_T = 0.0261032283758581;

    b2_T = -0.0982459868751664;
    c2_T = 0.20127470719351;
    d2_T = -0.173738749783189;

    b3_T = 0.156239779501715;
    c3_T = 0.0368592048084175;

    b4_T = -0.0478709034281346;

    CT = @(V,n,alpha,eps_eng)  a1_T + b1_T*J(V,n) + c1_T*J(V,n).^2 + d1_T*J(V,n).^3 + e1_T*J(V,n).^4 + ...
        phi(alpha,eps_eng).*(b2_T*J(V,n) + c2_T*J(V,n).^2 + d2_T*J(V,n).^3) +  phi(alpha,eps_eng).^2.*(b3_T*J(V,n) + c3_T*J(V,n).^2)...
        +  phi(alpha,eps_eng).^3.*b4_T*J(V,n);

    T = @(V,n,alpha,eps_eng) CT(V,n,alpha,eps_eng).*N_eng*rho*n.^2*D^4;

%X(1) = V
%X(2) = n
%X(3) = alpha
%X(4) = CD
%X(5) = CL
%X(6) = eps_eng

Vsol_fmincon = Xsol_fmincon(1);
nsol_fmincon = Xsol_fmincon(2);
alphasol_fmincon = Xsol_fmincon(3);
CLsol_fmincon = Xsol_fmincon(5);
CDsol_fmincon = Xsol_fmincon(4);
epsilonsol_fmincon = Xsol_fmincon(6);
% Constraint 1 : Hor forces

epscons1 = zeros(1,Nmat);
ncons1 = nlin;

%func1 = @(V) [ T(V,ncons1(jj),alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*V^2*S_ref*CD(alphasol_fmincon); 
   % T(X(1),X(2),alphasol_fmincon,epsilonsol_fmincon).*sin(phi(alphasol_fmincon,epsilonsol_fmincon))+0.5*rho*X(1).^2*S_ref*CL(alphasol_fmincon)-W];

for jj = 1:Nmat
    if jj == 1
        x0 = epsilonsol_fmincon;
    else
        x0 = epsilonsol_fmincon;
    end
epscons1( jj  )= fsolve(@(eps_eng) T(Vsol_fmincon,ncons1( jj ),alphasol_fmincon,eps_eng).*cos(phi(alphasol_fmincon,eps_eng))- 0.5*rho*Vsol_fmincon.^2*S_ref*CD(alphasol_fmincon), x0);
end
% disp("Valor falso?")
% T(24.20,14.54,alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*24.20.^2*S_ref*CD(alphasol_fmincon)
% disp("Valor verdadero?")
% T(24.2251,25.2146,alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*Vsol_fmincon.^2*S_ref*CD(alphasol_fmincon)
% Constraint 2: Ver Forces
epscons2 = zeros(1,Nmat);
ncons2 = nlin;


for jj = 1:Nmat 
    if jj == 1
        x0 = epsilonsol_fmincon;
    else
        x0 = epscons2(jj-1);
    end

   
epscons2( jj  )= fsolve(@(eps_eng) T(Vsol_fmincon,ncons2( jj   ),alphasol_fmincon,eps_eng).*sin(phi(alphasol_fmincon,eps_eng))+0.5*rho*Vsol_fmincon.^2*S_ref*CL(alphasol_fmincon)-W,x0);
end

% for jj = Nmat/2:-1:1
%     if jj == Nmat/2
%         x0 = epsilonsol_fmincon;
%     else
%         x0 = epscons2( jj+1 );
%     end
% epscons2(  Nmat/2 -jj +1 )= fsolve(@(eps_eng) T(Vsol_fmincon,ncons2( Nmat/2 -jj +1 ),alphasol_fmincon,eps_eng).*sin(phi(alphasol_fmincon,eps_eng))+0.5*rho*Vsol_fmincon.^2*S_ref*CL(alphasol_fmincon)-W,x0);
% end





epsplot2 = [epscons1; epscons2];
nplot2 = [ncons1; ncons2];


end

function [alphaplot3,Vplot3] = calculate_restrictions3 (Nmat, alphalin,Vlin,Xsol_fmincon)
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
D     = 0.8128; % Propeller Diameter [m]
RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;

%% Weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.5008; % Reference surface [m^2]
rho   = 1.2133;

%% Thrust


%% Propulsive Model

J   = @(V,n) V./(n*D);
phi = @(alpha, eps_eng) alpha + eps_eng;

 a1_T = 0.0735880531010883;
    b1_T = -0.0311758018412727;
    c1_T = -0.249744726429543;
    d1_T = 0.143084420694372;
    e1_T = 0.0261032283758581;

    b2_T = -0.0982459868751664;
    c2_T = 0.20127470719351;
    d2_T = -0.173738749783189;

    b3_T = 0.156239779501715;
    c3_T = 0.0368592048084175;

    b4_T = -0.0478709034281346;

    CT = @(V,n,alpha,eps_eng)  a1_T + b1_T*J(V,n) + c1_T*J(V,n).^2 + d1_T*J(V,n).^3 + e1_T*J(V,n).^4 + ...
        phi(alpha,eps_eng).*(b2_T*J(V,n) + c2_T*J(V,n).^2 + d2_T*J(V,n).^3) +  phi(alpha,eps_eng).^2.*(b3_T*J(V,n) + c3_T*J(V,n).^2)...
        +  phi(alpha,eps_eng).^3.*b4_T*J(V,n);

    T = @(V,n,alpha,eps_eng) CT(V,n,alpha,eps_eng).*N_eng*rho*n.^2*D^4;

%X(1) = V
%X(2) = n
%X(3) = alpha
%X(4) = CD
%X(5) = CL
%X(6) = eps_eng

Vsol_fmincon = Xsol_fmincon(1);
nsol_fmincon = Xsol_fmincon(2);
alphasol_fmincon = Xsol_fmincon(3);
CLsol_fmincon = Xsol_fmincon(5);
CDsol_fmincon = Xsol_fmincon(4);
epsilonsol_fmincon = Xsol_fmincon(6);
% Constraint 1 : Hor forces

alphacons1 = zeros(1,Nmat);
Vcons1 = Vlin;

%func1 = @(V) [ T(V,ncons1(jj),alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*V^2*S_ref*CD(alphasol_fmincon); 
   % T(X(1),X(2),alphasol_fmincon,epsilonsol_fmincon).*sin(phi(alphasol_fmincon,epsilonsol_fmincon))+0.5*rho*X(1).^2*S_ref*CL(alphasol_fmincon)-W];

for jj = 1:Nmat
    if jj == 1
        x0 = alphasol_fmincon;
    else
        x0 = alphacons1(jj-1);
    end
alphacons1( jj  )= fsolve(@(alpha) T(Vcons1(jj),nsol_fmincon,alpha,epsilonsol_fmincon).*cos(phi(alpha,epsilonsol_fmincon))- 0.5*rho*Vcons1(jj).^2*S_ref*CD(alpha), x0);
end
% disp("Valor falso?")
% T(24.20,14.54,alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*24.20.^2*S_ref*CD(alphasol_fmincon)
% disp("Valor verdadero?")
% T(24.2251,25.2146,alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*Vsol_fmincon.^2*S_ref*CD(alphasol_fmincon)
% Constraint 2: Ver Forces
alphacons2 = zeros(1,Nmat);
Vcons2 = Vlin;



for jj = 1:Nmat
    if jj == 1
        x0 = alphasol_fmincon;
    else
        x0 = alphacons2(jj-1);
    end

   
alphacons2( jj  )= fsolve(@(alpha) T(Vcons2(jj),nsol_fmincon,alpha,epsilonsol_fmincon).*sin(phi(alpha,epsilonsol_fmincon))+0.5*rho*Vcons2(jj).^2*S_ref*CL(alpha)-W,x0);
end

% for jj = Nmat/2:-1:1
%     if jj == Nmat/2
%         x0 = epsilonsol_fmincon;
%     else
%         x0 = epscons2( jj+1 );
%     end
% epscons2(  Nmat/2 -jj +1 )= fsolve(@(eps_eng) T(Vsol_fmincon,ncons2( Nmat/2 -jj +1 ),alphasol_fmincon,eps_eng).*sin(phi(alphasol_fmincon,eps_eng))+0.5*rho*Vsol_fmincon.^2*S_ref*CL(alphasol_fmincon)-W,x0);
% end





Vplot3 = [Vcons1; Vcons2];
alphaplot3 = [alphacons1; alphacons2];


end

function [alphaplot4,epsplot4] = calculate_restrictions4 (Nmat, alphalin,epsilonlin,Xsol_fmincon)
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
D     = 0.8128; % Propeller Diameter [m]
RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;

%% Weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.5008; % Reference surface [m^2]
rho   = 1.2133;

%% Thrust


%% Propulsive Model

J   = @(V,n) V./(n*D);
phi = @(alpha, eps_eng) alpha + eps_eng;

 a1_T = 0.0735880531010883;
    b1_T = -0.0311758018412727;
    c1_T = -0.249744726429543;
    d1_T = 0.143084420694372;
    e1_T = 0.0261032283758581;

    b2_T = -0.0982459868751664;
    c2_T = 0.20127470719351;
    d2_T = -0.173738749783189;

    b3_T = 0.156239779501715;
    c3_T = 0.0368592048084175;

    b4_T = -0.0478709034281346;

    CT = @(V,n,alpha,eps_eng)  a1_T + b1_T*J(V,n) + c1_T*J(V,n).^2 + d1_T*J(V,n).^3 + e1_T*J(V,n).^4 + ...
        phi(alpha,eps_eng).*(b2_T*J(V,n) + c2_T*J(V,n).^2 + d2_T*J(V,n).^3) +  phi(alpha,eps_eng).^2.*(b3_T*J(V,n) + c3_T*J(V,n).^2)...
        +  phi(alpha,eps_eng).^3.*b4_T*J(V,n);

    T = @(V,n,alpha,eps_eng) CT(V,n,alpha,eps_eng).*N_eng*rho*n.^2*D^4;

%X(1) = V
%X(2) = n
%X(3) = alpha
%X(4) = CD
%X(5) = CL
%X(6) = eps_eng

Vsol_fmincon = Xsol_fmincon(1);
nsol_fmincon = Xsol_fmincon(2);
alphasol_fmincon = Xsol_fmincon(3);
CLsol_fmincon = Xsol_fmincon(5);
CDsol_fmincon = Xsol_fmincon(4);
epsilonsol_fmincon = Xsol_fmincon(6);
% Constraint 1 : Hor forces

alphacons1 = zeros(1,Nmat);
epscons1 = epsilonlin;

%func1 = @(V) [ T(V,ncons1(jj),alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*V^2*S_ref*CD(alphasol_fmincon); 
   % T(X(1),X(2),alphasol_fmincon,epsilonsol_fmincon).*sin(phi(alphasol_fmincon,epsilonsol_fmincon))+0.5*rho*X(1).^2*S_ref*CL(alphasol_fmincon)-W];

for jj = 1:Nmat
    if jj == 1
        x0 = alphasol_fmincon;
    else
        x0 = alphacons1(jj-1);
    end
alphacons1( jj  )= fsolve(@(alpha) T(Vsol_fmincon,nsol_fmincon,alpha,epscons1(jj)).*cos(phi(alpha,epscons1(jj)))- 0.5*rho*Vsol_fmincon.^2*S_ref*CD(alpha), x0);
end
% disp("Valor falso?")
% T(24.20,14.54,alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*24.20.^2*S_ref*CD(alphasol_fmincon)
% disp("Valor verdadero?")
% T(24.2251,25.2146,alphasol_fmincon,epsilonsol_fmincon).*cos(phi(alphasol_fmincon,epsilonsol_fmincon))- 0.5*rho*Vsol_fmincon.^2*S_ref*CD(alphasol_fmincon)
% Constraint 2: Ver Forces
alphacons2 = zeros(1,Nmat);
epscons2 = epsilonlin;



for jj = 1:Nmat
    if jj == 1
        x0 = alphasol_fmincon;
    else
        x0 = alphacons2(jj-1);
    end

   
alphacons2( jj  )= fsolve(@(alpha) T(Vsol_fmincon,nsol_fmincon,alpha,epscons2(jj)).*sin(phi(alpha,epscons2(jj)))+0.5*rho*Vsol_fmincon.^2*S_ref*CL(alpha)-W,x0);
end

% for jj = Nmat/2:-1:1
%     if jj == Nmat/2
%         x0 = epsilonsol_fmincon;
%     else
%         x0 = epscons2( jj+1 );
%     end
% epscons2(  Nmat/2 -jj +1 )= fsolve(@(eps_eng) T(Vsol_fmincon,ncons2( Nmat/2 -jj +1 ),alphasol_fmincon,eps_eng).*sin(phi(alphasol_fmincon,eps_eng))+0.5*rho*Vsol_fmincon.^2*S_ref*CL(alphasol_fmincon)-W,x0);
% end





epsplot4 = [epscons1; epscons2];
alphaplot4 = [alphacons1; alphacons2];


end