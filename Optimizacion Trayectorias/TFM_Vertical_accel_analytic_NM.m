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


CT0 =  0.0735880531010883;
CT1 = -0.0311758018412727;
CT2 = -0.249744726429543;  


CP0 = 0.0261518307541734;
CP1 = 0.0473735972985378;
CP2 = -0.16267474946046;
CP3 = 0.0247028469343899;

%% Loads and weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S     = 0.5008; % Reference surface [m^2]
rho   = 1.2133;



%% Electric
tau = 0.2;
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency
mbatt = 3; %[kg]
e0    = 720e3; %[J/kg]
E     = e0*mbatt;  %[J] Total energy of the battery packs
N_eng = 2; % Number of engines
D     = 0.8128; % Propeller Diameter [m]

RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;

%% Thrust
CT = @(V,n) CT2*V^2/(n^2*D^2)+CT1*V/(n*D)+CT0;

T = @(V,n) CT(V,n)*N_eng*rho*n^2*D^4;

%% Power
CP = @(V,n) CP3*V^3/(n^3*D^3)+CP2*V^2/(n^2*D^2)+CP1*V/(n*D)+CP0;

P = @(V,n) CP(V,n)*N_eng*rho*n^3*D^5;


%% Set deltaH (altitude difference) 

deltaH = 20; %[m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First Method: Pseudo-Analytical with integrals:


% First, we define the constants:

a_int = rho*(N_eng*D^2*CT2-1/2*S*CD(alpha_V))/mTOW;
b_int = @(n) rho*(N_eng*n*D^3*CT1)/mTOW;
c_int = @(n) rho*(N_eng*n^2*D^4*CT0)/mTOW-g;



% Define vector for rps and deltaH:
Nlin = 20;
Hlin = 150;
n_v      = linspace(50,rps_max,Nlin);
% Start interpolation from known segment (IDEA)

deltaH_v = linspace(0,deltaH,Hlin);

%% Preallocate for speed
V_sol = zeros(Nlin,1);

H_v   = zeros(Hlin,1);
t_v   = zeros(Hlin,Nlin);
P_v   = zeros(Hlin,Nlin);
T_v   = zeros(Hlin,Nlin);
nflag = zeros(Hlin,1);
Vflag = zeros(Hlin,1);
Hflag = zeros(Hlin,1);
tflag = zeros(Hlin,1);
iiflag = zeros(Hlin,1);
jjflag = zeros(Hlin,1);



V_sol(1) = 0;
H_v  (1) = 0;
t_v  (1) = 0;
options = optimset('TolX',1e-5,'MaxIter',2000);

%% Start loop
for ii = 1:Nlin
    % Define revolutions at which we are flying
    n_input = n_v(ii);
    
    % Start climb from H = 0
    V00 = 0;
    V_0 = 0.5;
    
    % This flag is for calculating the point where the Power limitation is
    % exceeded, if it exists.
    exceedflag = 0;
    kk = 1;
    
    % To calculate Power at first instant
    P_v(1,ii) = P(0,n_input);
    T_v(1,ii) = T(0,n_input);
    
    for jj = 2:Hlin
        
        h_step = deltaH_v(jj)-deltaH_v(jj-1);
        
        % We calculate the analytical integral:
        % V_ev is the speed we are evaluating the function.
        fV_ev = @(V_ev,n) log(abs(a_int*V_ev^2+b_int(n)*V_ev+c_int(n)))/(2*a_int)-(b_int(n)*atan((2*a_int*V_ev+b_int(n))/sqrt(4*a_int*c_int(n)-b_int(n)^2)))...
            /(a_int*sqrt(4*a_int*c_int(n)-b_int(n)^2)) ;
        
        
        
        % Solve VF for a given RPM
        
        
        V_sol(jj,ii) = fsolve(@(V_ev) fV_ev(V_ev,n_input)- fV_ev(V00,n_input) - h_step,V_0,options);
        
        P_v (jj,ii) = P(V_sol(jj,ii),n_input);
        T_v (jj,ii) = T(V_sol(jj,ii),n_input);
        
        
        % For the next step, we start from a higher speed
        V00 = V_sol(jj,ii);
        % Update initial iterant
        V_0 = V_sol(jj,ii);
        H_v (jj) = h_step + H_v(jj-1);
        t_v (jj,ii) = h_step/V_sol(jj,ii) + t_v(jj-1,ii);
        
        %% We account if the mission is feasible in terms of Thrust
        if T_v(jj,ii) > Tmax && exceedflag == 0
            nflag(kk) = n_input;
            Vflag(kk) = V_sol(jj,ii);
            Hflag(kk) = H_v (jj);
            tflag(kk) = t_v (kk);
            iiflag(kk) = ii;
            jjflag(kk) = jj;
            exceedflag = 1;
            kk = kk+1;
            warning('THRUST has been Exceeded')
        end
        
    end
    
end



%% Calculate Energy for each iteration

VF_sol = V_sol(end,:)';

% Define the constants of the Power function

ap_int = rho*N_eng*CP3*D^2/eta_m/(1-tau);
bp_int = @(n) rho*N_eng*CP2*n*D^3/eta_m/(1-tau);
cp_int = @(n) rho*N_eng*CP1*n^2*D^4/eta_m/(1-tau);
dp_int = @(n) rho*N_eng*CP0*n^3*D^5/eta_m/(1-tau);


f_Energy = @(V_E,n)-((-a_int^2*cp_int(n)+a_int*b_int(n)*bp_int(n)+a_int*c_int(n)*ap_int-b_int(n)^2*ap_int)*log(abs(a_int*V_E^2+b_int(n)*V_E+c_int(n))))/...
    (2*a_int^3)-((-a_int^2*dp_int(n)-(b_int(n)*(-a_int^2*cp_int(n)+a_int*b_int(n)*bp_int(n)+a_int*c_int(n)*ap_int-b_int(n)^2*ap_int))/(2*a_int)+c_int(n)*(a_int*bp_int(n)-b_int(n)*ap_int))...
    *atan((2*a_int*V_E+b_int(n))/(2*sqrt(a_int)*sqrt(c_int(n)-b_int(n)^2/(4*a_int)))))/(a_int^(5/2)*sqrt(c_int(n)-b_int(n)^2/(4*a_int)))...
    +(ap_int*V_E^2)/(2*a_int)+(bp_int(n)/a_int-(b_int(n)*ap_int)/a_int^2)*V_E;

E_v = zeros(Nlin,1);
for ii = 1:Nlin
    
    E_v (ii) =   f_Energy(VF_sol(ii),n_v(ii))-f_Energy(0,n_v(ii));
    
    
end

[minEv,indexminE] = min(E_v);

%% Calculate optimum rps solution

% RPS that give minimum Energy consumption
n_min = n_v(indexminE);


%% Define matrices for plot
Hmat = repmat(H_v',Nlin,1)';
nmat = repmat(n_v,Hlin,1);
tmat = repmat(t_v',Nlin,1)';

vect_cc_n = linspace(min(min(nmat)),max(max(nmat)),Nlin)';


figure(1)
[n_c,h_nmat_c] = contourf([Hmat(:,1:indexminE-1) Hmat(:,indexminE+1:Nlin)],...
    [V_sol(:,1:indexminE-1) V_sol(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'k','LineWidth',1.5);
clabel(n_c,h_nmat_c)
hold on
plot (Hmat(:,indexminE),V_sol(:,indexminE),'r','LineWidth',1.5);
plot (Hflag,Vflag,'sm');
colormap(flipud(colormap('white')))
grid on
Title1 = strcat(['Vertical speed [m/s] vs Height [m] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
title(Title1)
xlabel('Climb height [m] ')
ylabel('V [m/s]')
caxis([min(min(nmat)) max(max(nmat))]); %%%%%%%%%%%%%%%



figure(2)

plot(n_v,E_v,'k','LineWidth',1.8)
grid on
hold on
plot(n_min,minEv,'or','LineWidth',1.2)
xline(rps_max,'--b','LineWidth',1.5)
xlabel('Engine revolutions [rps] ')
ylabel('Energy consumption [J]')
Title2 = strcat('Energy consumption vs Engine rev. [rps]');
title(Title2)
legend('E(n)',['Opt Energy consumption = ' num2str(minEv/1000) ' kJ.'],'Max RPS','Location','best')

figure(3)

[n_c2,t_nmat_c] = contourf([tmat(:,1:indexminE-1) tmat(:,indexminE+1:Nlin)],...
    [V_sol(:,1:indexminE-1) V_sol(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'k','LineWidth',1.5);
clabel(n_c2,t_nmat_c)
colormap(flipud(colormap('white')))
hold on
plot (tmat(:,indexminE),V_sol(:,indexminE),'r','LineWidth',1.5);
plot (tflag,Vflag,'sm','LineWidth',1.5);
grid on
Title3 = strcat(['Vertical speed [m/s] vs time [s] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
title(Title3)
xlabel('Climb time [s] ')
ylabel('V [m/s]')
caxis([min(min(nmat)) max(max(nmat))]); %%%%%%%%%%%%%%%


figure(4)

[n_c3,ht_nmat_c] = contourf([tmat(:,1:indexminE-1) tmat(:,indexminE+1:Nlin)],...
    [Hmat(:,1:indexminE-1) Hmat(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'k','LineWidth',1.5);
clabel(n_c3,ht_nmat_c)
colormap(flipud(colormap('gray')))
hold on
plot (tmat(:,indexminE),Hmat(:,indexminE),'r','LineWidth',1.5);
plot (tflag,Hflag,'sm','LineWidth',1.5);

grid on
Title3 = strcat(['Climb height [m] vs time [s] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
title(Title3)
xlabel('Climb time [s] ')
ylabel('Climb height [m]')
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
caxis([min(min(nmat)) 1.2*max(max(nmat))]); %%%%%%%%%%%%%%%

figure(5)
[n_c,P_nmat_c] = contourf([Hmat(:,1:indexminE-1) Hmat(:,indexminE+1:Nlin)],...
    [P_v(:,1:indexminE-1) P_v(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'k','LineWidth',1.5);
clabel(n_c,P_nmat_c)
hold on
plot (Hmat(:,indexminE),P_v(:,indexminE),'r','LineWidth',1.5);
%yline(Pmax,'-.b','LineWidth',2);
colormap(flipud(colormap('white')))
grid on
Title1 = strcat(['Power consumption [W] vs Height [m] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
title(Title1)
xlabel('Climb height [m] ')
ylabel('Power [W]')
caxis([min(min(nmat)) max(max(nmat))]); %%%%%%%%%%%%%%%

