%% Range maximization for Cruise performance
%% Now, you can chose the aerodynamic model!


%% Code for solving Max Range for a given Energy
clear all
close all

%% Two methods, graphic and using fmincon
% Analytical methodologies are not available because we have an order 9
% polynomial on CL and CD


%% Aerodynamic model

% Set AeroModel = 1 for the full coefficient, alpha dependant model.
% Set AeroModel = 0 for the drag polar CD0-k model.

%%%%%%%%%%%%%%%%%%%%%%%
AeroModel = 1; % <---------------------------
%%%%%%%%%%%%%%%%%%%%%%%


if AeroModel == 1
    CD = @(alpha) 0.03356*(alpha/52.18).^9-0.02777*(alpha/52.18).^8-0.2525*(alpha/52.18).^7+0.3034*(alpha/52.18).^6+...
        0.6575*(alpha/52.18).^5-1.635*(alpha/52.18).^4-0.6442*(alpha/52.18).^3+3.082*(alpha/52.18).^2+0.01637*(alpha/52.18)+0.08034;
    
    
    CL = @(alpha) 0.02081*(alpha/52.18).^8-0.002541*(alpha/52.18).^7-0.1681*(alpha/52.18).^6+0.2368*(alpha/52.18).^5+...
        0.491*(alpha/52.18).^4-1.652*(alpha/52.18).^3-0.533*(alpha/52.18).^2+3.047*(alpha/52.18)+0.582;
elseif  AeroModel == 0
    
    CL_max_w1_CR  = 1.5008;
    Safety_F      = 1.2;
    CL_max_w1_CR_ope = CL_max_w1_CR/Safety_F;
    
    CD0 = 0.035682195723975;
    CD2 = 0.054209627025009;
    
end


%% Propulsive

CT0 =  0.089050757500461;
CT1 = -0.026659496766654;
CT2 = -0.162620845549081;


CP0 = 0.034214580122684;
CP1 = -0.002523708791417;
CP2 = 0.116121898742278;
CP3 = -0.248807360672063;

%% Loads and Weight
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Firstly, solve in a graphic manner

%% Power
CP = @(V,n) CP3*V^3/(n^3*D^3)+CP2*V^2/(n^2*D^2)+CP1*V/(n*D)+CP0;

P = @(V,n) CP(V,n)*N_eng*rho*n^3*D^5;


%% Range

xf = @(V,n)  eta_m*E*V/((1-tau)*P(V,n))*1e-3;



%% Generate numerical Matrix
Nmat = 100;
Vlin = linspace(18,40,Nmat);
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

if AeroModel ==1
    %% Calculate restriction L = W 
    CLcon = @(V) W./(1/2*rho*V.^2*S_ref); 
    
    ncon = @(V,alpha) -(1/2)*(D*N_eng*CT1-sqrt(-4*CT0*CT2*D^2*N_eng^2+CT1^2*D^2*N_eng^2+2*CD(alpha)*CT0*N_eng*S_ref))*V/(CT0*N_eng*D^2);
    % Preallocate vectors for speed
    nconlin = zeros(Nmat,1);
    alphaconlin = zeros(Nmat,1);
    
    alpha0 = 6; % First iteration
    for kk = 1:Nmat
        
        % Solve the angle that makes the CL that the speed V imposes
        
        alphaconlin(kk) = fzero(@(alpha) -CLcon(Vlin(kk))+CL(alpha),alpha0);
        
        % With the calculated angle, calculate D, then solve revolutions n
        nconlin (kk) = ncon(Vlin(kk),alphaconlin(kk));
        
        
    end
    
    %% Vstall calculations (for plot)
    alphavec = linspace(-90,90,100);
    CLvec = CL(alphavec);
    CLmax = max(CLvec);
    
    Vstall = sqrt(2*W/(rho*S_ref*CLmax));
    
elseif AeroModel ==0
    ncon = @(V) -(1/2)*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+...
        8*CT0*N_eng*S_ref*W^2*CD2))/(CT0*N_eng*S_ref*rho*V*D^2);
    % Preallocate matrix for speed
    nconlin = zeros(Nmat,1);
    for kk = 1:Nmat
        nconlin (kk) = ncon(Vlin(kk));
        
    end
    %% Calculate stall speed
    
    Vstall = sqrt(2*W/(rho*S_ref*CL_max_w1_CR_ope));
end




%% Eliminate zeros from indeterminations (points where CP = 0)
xfmat  = xfmat';
threshold = 2e2; % Assume Range < 200 km 
xfmat(xfmat > threshold) = 0;
xfmat(xfmat < 0) = 0;




%% Now, solve using fmincon




if AeroModel ==1
 %X(1) = V
%X(2) = n
%X(3) = alpha
%X(4) = CD
%X(5) = CL




%% Power
CP_fmin = @(X) CP3*X(1)^3/(X(2)^3*D^3)+CP2*X(1)^2/(X(2)^2*D^2)+CP1*X(1)/(X(2)*D)+CP0;

P_fmin = @(X) CP_fmin(X)*N_eng*rho*X(2)^3*D^5;
   
  %% Range

% - because fmincon solves for the minimum of a function!!
xf_fmin = @(X) -eta_m*E*X(1)/((1-tau)*P_fmin(X))*1e-3;
nonlcon = @constrains_TD_1;  
    
x0 = [21.08 50.9 4 0.08 0.8];
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
Xsol = fmincon(xf_fmin,x0,A,b,Aeq,beq,lb,ub,nonlcon);

nsol_fmincon = Xsol(2);
Vsol_fmincon = Xsol(1);


elseif AeroModel ==0
    
%% Power
CP_f = @(X) CP3*X(1)^3/(X(2)^3*D^3)+CP2*X(1)^2/(X(2)^2*D^2)+CP1*X(1)/(X(2)*D)+CP0;

P_f = @(X) CP_f(X)*N_eng*rho*X(2)^3*D^5;

%% Range

% - because fmincon solves for the minimum of a function!!
xf_fmin = @(X) -eta_m*E*X(1)/((1-tau)*P_f(X))*1e-3;


nonlcon = @constrains_TD_0;


x0_f = [24.08 55.9];
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];   
    
Xsol = fmincon(xf_fmin,x0_f,A,b,Aeq,beq,lb,ub,nonlcon);
nsol_fmincon = Xsol(2);
Vsol_fmincon = Xsol(1); 
    
    
    
end

%% This maximum considers the stall limitation
xf_max_fmin       = -xf_fmin (Xsol);



%%%%%%
% Now everything is solved!! We only need to plot the results
%%%%%%


%% Plot contours
N_contour_lines = 10; % Number of contour lines
vect_cc_xf = linspace(min(min(xfmat)),max(max(xfmat)),N_contour_lines);

%% Plot results
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
         xline(rps_max,':b','LineWidth',2)
         plot(nsol_fmincon,Vsol_fmincon,'om','LineWidth',2.5)
         legend('x_f (V,n)','T = D constrain','Stall Speed','Max RPS','Optimal solution','Location','northwest')



%% Auxiliary function for fmincon

function [c,ceq] = constrains_TD_0(X)
%% Aerodynamic model
CL_max_w1_CR  = 1.5008;
Safety_F      = 1.2;
CL_max_w1_CR_ope = CL_max_w1_CR/Safety_F;



N_eng = 2; % Number of engines
D     = 0.7112; % Propeller Diameter [m]
RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter



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

% [Stall speed condition; Max RPS condition]s
c  = [+sqrt(2*W/(rho*S_ref*CL_max_w1_CR_ope))-X(1);-rps_max+X(2)]; 

% T = D condition
ceq = N_eng*rho*(CT2*X(1)^2/(X(2)^2*D^2)+CT1*X(1)/(X(2)*D)+CT0)*X(2)^2*D^4 ... 
    -(1/2)*rho*X(1)^2*S_ref*(4*W^2*CD2/(rho^2*X(1)^4*S_ref^2)+CD0);
end


function [c,ceq] = constrains_TD_1(X)


%% Propulsive model
N_eng = 2; % Number of engines
D     = 0.7112; % Propeller Diameter [m]
RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter


CT0 =  0.089050757500461;
CT1 = -0.026659496766654;
CT2 = -0.162620845549081; 

%% Loads and weights

mTOW = 21.39; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.430279179101; % Reference surface [m^2]
rho   = 1.225;

%% Vstall calculations

CL = @(alpha) 0.02081*(alpha/52.18).^8-0.002541*(alpha/52.18).^7-0.1681*(alpha/52.18).^6+0.2368*(alpha/52.18).^5+...
    0.491*(alpha/52.18).^4-1.652*(alpha/52.18).^3-0.533*(alpha/52.18).^2+3.047*(alpha/52.18)+0.582;


alphavec = linspace(-90,90,100);
CLvec = CL(alphavec);
CLmax = max(CLvec);

Vstall = sqrt(2*W/(rho*S_ref*CLmax));

%X(1) = V
%X(2) = n
%X(3) = alpha
%X(4) = CD
%X(5) = CL


% Stall speed condition
c  = [Vstall-X(1);-rps_max+X(2)];


% T = D condition
ceq =  [N_eng*rho*(CT2*X(1)^2/(X(2)^2*D^2)+CT1*X(1)/(X(2)*D)+CT0)*X(2)^2*D^4-(1/2)*rho*X(1)^2*S_ref*X(4);
    sqrt(2*W/(rho*S_ref*X(5)))-X(1);
    -X(4)+0.03356*(X(3)/52.18)^9-0.02777*(X(3)/52.18)^8-0.2525*(X(3)/52.18)^7+0.3034*(X(3)/52.18)^6+0.6575*(X(3)/52.18)^5-1.635*(X(3)/52.18)^4-0.6442*(X(3)/52.18)^3+3.082*(X(3)/52.18)^2+0.01637*(X(3)/52.18)+0.08034;
    -X(5)+0.02081*(X(3)/52.18)^8-0.002541*(X(3)/52.18)^7-0.1681*(X(3)/52.18)^6+0.2368*(X(3)/52.18)^5+0.491*(X(3)/52.18)^4-1.652*(X(3)/52.18)^3-0.533*(X(3)/52.18)^2+3.047*(X(3)/52.18)+0.582];


end





