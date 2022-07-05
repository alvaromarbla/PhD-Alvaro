%% Code for solving Min Enegy consumption for a altitude difference for Vertical Take-off
clear all
close all




%% Aerodynamic model
alpha_V = -90*pi/180; % [º] Because of Vertical TO

p_CL_ac = [0.582;3.345731;-0.642635;-2.187085;0.713766;0.377985;-0.2946314;-0.004890;0.043976];


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


% Full model w.r.t angle of attack
CL = @(alpha) p_CL_ac(9)*alpha.^8+p_CL_ac(8)*alpha.^7 + p_CL_ac(7)*alpha.^6 + p_CL_ac(6)*alpha.^5 + p_CL_ac(5)*alpha.^4 + p_CL_ac(4)*alpha.^3 ...
    + p_CL_ac(3)*alpha.^2 + p_CL_ac(2)*alpha + p_CL_ac(1);

CD = @(alpha) p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + p_CD_ac(4)*alpha.^3 ...
    + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);
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

S = 0.430279179101; % Reference surface [m^2]
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

deltaH = 20; %[m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% We start with graphic method 


%% Power
CP = @(V,n) CP3*V^3/(n^3*D^3)+CP2*V^2/(n^2*D^2)+CP1*V/(n*D)+CP0;

P = @(V,n) CP(V,n)*N_eng*rho*n^3*D^5;

%% Energy

E_v = @(V,n)  P(V,n)*deltaH/(eta_m*V); % -> Solve dE/dV = 0?

%% Generate numerical Matrix
Nmat = 1000;
Vlin = linspace(0,30,Nmat);
nlin = linspace(0,110, Nmat);

% Preallocate matrix for speed
Ev_mat = zeros(Nmat,Nmat);
% Loop in speed
for  ii = 1: Nmat
    % Loop in revolutions
    for jj = 1:Nmat
        
        Ev_mat(ii,jj) = E_v(Vlin(jj),nlin(ii));
        
    end
end

%% Calculate restriction T = D + W

ncon = @(V,alpha) -(1/2)*(CT1*N_eng*V*rho*D-sqrt(-4*CT0*CT2*D^2*N_eng^2*V^2*rho^2+CT1^2*D^2*N_eng^2*V^2*rho^2+2*CD(alpha)*CT0*N_eng*S*V^2*rho^2+4*CT0*N_eng*W*rho))/(CT0*N_eng*rho*D^2)
% Preallocate matrix for speed
nconlin = zeros(Nmat,1);


for kk = 1:Nmat
    nconlin (kk) = ncon(Vlin(kk),alpha_V);

end


%% Eliminate zeros from indeterminations (points where CP = 0)
Ev_mat  = Ev_mat'; % To make revolutions stand on the x-axis and V on the y-axis


threshold = 2e4; % Assume Range < 200 km 
Ev_mat(Ev_mat > threshold) = 0;
Ev_mat(Ev_mat < 0) = 0;


 
%% Now, solve in an Analytical manner

% It is mandatory to evaluate CD as CD_alpha = 90º

CD_V = CD(alpha_V);


%% Function Definition
%Derivative of xf w.r.t V
fsol = @(Vv) -(3/8).*(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^2.*...
    (-8.*CP3.*Vv^3.*CT0^3.*N_eng^3.*rho^3.*D^3/(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^3 ...
+4.*CP2.*Vv^2.*CT0^2.*N_eng^2.*rho^2.*D^2/(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^2-...
2.*CP1.*Vv.*CT0.*N_eng.*rho.*D/(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))+CP0).*deltaH.*(CT1.*N_eng.*rho.*D-(1/2).*...
(-8.*CT0.*CT2.*D^2.*N_eng^2.*Vv.*rho^2+2.*CT1^2.*D^2.*N_eng^2.*Vv.*rho^2+4.*CD_V.*CT0.*N_eng.*S.*Vv.*rho^2)/sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))...
/(N_eng^2.*rho^2.*CT0^3.*D.*Vv)-(1/8).*(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^3 ...
.*(-24.*CP3.*Vv^2.*CT0^3.*N_eng^3.*rho^3.*D^3/(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^3 ...
+24.*CP3.*Vv^3.*CT0^3.*N_eng^3.*rho^3.*D^3.*(CT1.*N_eng.*rho.*D-(1/2).*(-8.*CT0.*CT2.*D^2.*N_eng^2.*Vv.*rho^2+2.*CT1^2.*D^2.*N_eng^2.*Vv.*rho^2+4.*CD_V.*CT0.*N_eng.*S.*Vv.*rho^2)/...
sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))/(CT1.*N_eng.*Vv.*rho.*D-...
sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^4+8.*CP2.*Vv.*CT0^2.*N_eng^2.*rho^2.*D^2/...
(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^2 ...
-8.*CP2.*Vv^2.*CT0^2.*N_eng^2.*rho^2.*D^2.*(CT1.*N_eng.*rho.*D-(1/2).*(-8.*CT0.*CT2.*D^2.*N_eng^2.*Vv.*rho^2+2.*CT1^2.*D^2.*N_eng^2.*Vv.*rho^2+4.*CD_V.*CT0.*N_eng.*S.*Vv.*rho^2)...
/sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))/(CT1.*N_eng.*Vv.*rho.*D-...
sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^3-2.*CP1.*CT0.*N_eng.*rho.*D/...
(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))+...
2.*CP1.*Vv.*CT0.*N_eng.*rho.*D.*(CT1.*N_eng.*rho.*D-(1/2).*(-8.*CT0.*CT2.*D^2.*N_eng^2.*Vv.*rho^2+2.*CT1^2.*D^2.*N_eng^2.*Vv.*rho^2+4.*CD_V.*CT0.*N_eng.*S.*Vv.*rho^2)/...
sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))/(CT1.*N_eng.*Vv.*rho.*D-...
sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^2).*deltaH/(N_eng^2.*rho^2.*CT0^3.*D.*Vv)...
+(1/8).*(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^3.*...
(-8.*CP3.*Vv^3.*CT0^3.*N_eng^3.*rho^3.*D^3/(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^3+...
4.*CP2.*Vv^2.*CT0^2.*N_eng^2.*rho^2.*D^2/(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^2-...
2.*CP1.*Vv.*CT0.*N_eng.*rho.*D/(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))+CP0).*deltaH/...
(N_eng^2.*rho^2.*CT0^3.*D.*Vv^2).*((CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^3.*...
(-8.*CP3.*Vv^3.*CT0^3.*N_eng^3.*rho^3.*D^3/(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^3+...
4.*CP2.*Vv^2.*CT0^2.*N_eng^2.*rho^2.*D^2/(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))^2-...
2.*CP1.*Vv.*CT0.*N_eng.*rho.*D/(CT1.*N_eng.*Vv.*rho.*D-sqrt(-4.*CT0.*CT2.*D^2.*N_eng^2.*Vv^2.*rho^2+CT1^2.*D^2.*N_eng^2.*Vv^2.*rho^2+2.*CD_V.*CT0.*N_eng.*S.*Vv^2.*rho^2+4.*CT0.*N_eng.*W.*rho))+CP0).*deltaH);




%% Set fsolve
options = optimset('TolX',1e-5,'MaxIter',2000);


X0 = 15;
[Vsol,err] = fsolve(fsol,X0,options);



% Calculate the revolutions from the T = D + W restriction
nsol_A = ncon(Vsol,alpha_V);


%% Finally, solve using fmincon

%% Power
CP_f = @(X) CP3*X(1)^3/(X(2)^3*D^3)+CP2*X(1)^2/(X(2)^2*D^2)+CP1*X(1)/(X(2)*D)+CP0;

P_f = @(X) CP_f(X)*N_eng*rho*X(2)^3*D^5;

%% Range

E_v_fmincon = @(X) P_f(X)*deltaH/(eta_m*X(1));  

nonlcon = @constrains_TD;


x0_f = [15 60];
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];


Xsol_fmincon = fmincon(E_v_fmincon,x0_f,A,b,Aeq,beq,lb,ub,nonlcon)
nsol_fmincon = Xsol_fmincon(2);
Vsol_fmincon = Xsol_fmincon(1);


%%  Minimum evaluation
E_min_fmin       = E_v_fmincon (Xsol_fmincon);



%%%%%%
% Now everything is solved!! We only need to plot the results
%%%%%%

%% Plot contours
N_contour_lines = 12; % Number of contour lines
vect_cc_E = linspace(min(min(Ev_mat)),max(max(Ev_mat)),N_contour_lines);


 figure(1)
[E_c,h_Emat_c] = contourf(nlin,Vlin,Ev_mat,vect_cc_E');
       clabel(E_c,h_Emat_c)
       colormap(flipud(colormap('gray')))
         grid on
         Title1 = strcat(['Energy consumption [J] vs. V [m/s] & Engine rev. [rps]. Min Energy is ', num2str(E_min_fmin), ' J']);
         title(Title1)
         xlabel('Engine revolutions [rps] ')
         ylabel('V [m/s]')
         sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
         sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
         sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
         sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
        caxis([min(min(Ev_mat)) max(max(Ev_mat))]); %%%%%%%%%%%%%%%
         hold on 
         plot(nconlin,Vlin,'--r','LineWidth',2) % plot constrain
         xline(rps_max,':b','LineWidth',2)
         plot(nsol_A,Vsol,'*g','LineWidth',2)
         plot(nsol_fmincon,Vsol_fmincon,'om','LineWidth',2.5)
         legend('E (V,n)','T = D + W constrain','Max RPS','Optimal solution Analytical','Optimal solution Fmincon','Location','northwest')

         battery_vo = E_min_fmin/e0



%% Auxiliary function for fmincon
function [c,ceq] = constrains_TD(X)
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

S = 0.430279179101; % Reference surface [m^2]
rho   = 1.225;

%% Drag model
alpha_V = -90; % [º] Because of Vertical TO

CD = @(alpha) 0.03356*(alpha/52.18).^9-0.02777*(alpha/52.18).^8-0.2525*(alpha/52.18).^7+0.3034*(alpha/52.18).^6+...
        0.6575*(alpha/52.18).^5-1.635*(alpha/52.18).^4-0.6442*(alpha/52.18).^3+3.082*(alpha/52.18).^2+0.01637*(alpha/52.18)+0.08034;

CD_V = CD(alpha_V);
%X(1) = V
%X(2) = n

% [Stall speed condition; Max RPS conditions]



c  = [-X(1);-X(2);-rps_max+X(2)]; 

% T = D condition
ceq = N_eng*rho*X(2)^2*D^4*(CT2*X(1)^2/(X(2)^2*D^2)+CT1*X(1)/(X(2)*D)+CT0)-(1/2)*rho*X(1)^2*S*CD_V-W;
end
