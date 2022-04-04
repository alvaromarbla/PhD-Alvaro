%% Code for solving Max Range for a given Energy
clear all
close all




%% Aerodynamic model
CL_max_w1_CR  = 1.5008;
Safety_F      = 1.2;
CL_max_w1_CR_ope = CL_max_w1_CR/Safety_F;

CD0 = 0.035682195723975;
CD2 = 0.054209627025009;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% We start with graphic method 


%% Power
CP = @(V,n) CP3*V^3/(n^3*D^3)+CP2*V^2/(n^2*D^2)+CP1*V/(n*D)+CP0;

P = @(V,n) CP(V,n)*N_eng*rho*n^3*D^5;

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

ncon = @(V) -(1/2)*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+...
    8*CT0*N_eng*S_ref*W^2*CD2))/(CT0*N_eng*S_ref*rho*V*D^2);
% Preallocate matrix for speed
nconlin = zeros(Nmat,1);
for kk = 1:Nmat
    nconlin (kk) = ncon(Vlin(kk));

end

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
[Vsol,err] = fsolve(fsol,X0,options);

ncon_A = @(V) -(1/2)*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+...
    8*CT0*N_eng*S_ref*W^2*CD2))/(CT0*N_eng*S_ref*rho*V*D^2);

nsol_A = ncon_A(Vsol);


%% Finally, solve using fmincon

%% Power
CP_f = @(X) CP3*X(1)^3/(X(2)^3*D^3)+CP2*X(1)^2/(X(2)^2*D^2)+CP1*X(1)/(X(2)*D)+CP0;

P_f = @(X) CP_f(X)*N_eng*rho*X(2)^3*D^5;

%% Range

% - because fmincon solves for the minimum of a function!!
xf_fmincon = @(X) -eta_m*E*X(1)/((1-tau)*P_f(X))*1e-3;
nonlcon = @constrains_TD;


x0_f = [24.08 55.9];
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

%% Plot contours
N_contour_lines = 10; % Number of contour lines
vect_cc_xf = linspace(min(min(xfmat)),max(max(xfmat)),N_contour_lines);


 figure(1)
[xf_c,h_xfmat_c] = contourf(nlin,Vlin,xfmat,vect_cc_xf');
       clabel(xf_c,h_xfmat_c)
       colormap(flipud(colormap('gray')))
         grid on
         Title1 = strcat(['Range [km] vs. V [m/s] & Engine rev. [rps]. Max Range is ', num2str(xf_max_fmin), ' km']);
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
         yline(Vstall,'-.b','LineWidth',2)
         plot(nsol_A,Vsol,'*g','LineWidth',2)
         plot(nsol_fmincon,Vsol_fmincon,'om','LineWidth',2.5)
         legend('x_f (V,n)','T = D constrain','Stall Speed','Optimal solution no-stall','Optimal solution stall','Location','northwest')

% figure(2)
% 
% mesh(Vlin,nlin,xfmat)



%% Auxiliary function for fmincon
function [c,ceq] = constrains_TD(X)
%% Aerodynamic model
CL_max_w1_CR  = 1.5008;
Safety_F      = 1.2;
CL_max_w1_CR_ope = CL_max_w1_CR/Safety_F;



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

%X(1) = V
%X(2) = n

% Stall speed condition
c  = sqrt(2*W/(rho*S_ref*CL_max_w1_CR_ope))-X(1); 

% T = D condition
ceq = N_eng*rho*(CT2*X(1)^2/(X(2)^2*D^2)+CT1*X(1)/(X(2)*D)+CT0)*X(2)^2*D^4 ... 
    -(1/2)*rho*X(1)^2*S_ref*(4*W^2*CD2/(rho^2*X(1)^4*S_ref^2)+CD0);
end
