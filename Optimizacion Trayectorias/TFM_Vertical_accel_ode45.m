%% Code for solving Min Enegy consumption for a altitude difference for Vertical Take-off

%% Use of ode45 commands
clear all
%close all
%% Set deltaH (altitude difference)
global deltaH
deltaH = 100; %[m]
global D
D     = 0.7112; % Propeller Diameter [m]
global N_eng;
N_eng = 2; % Number of engines
global rho
rho   = 1.225; % Air density [kg/m^3]

%% Define maximum rps


RPMMAX_APC = 150000; % Max RPM of AXI motors
D_inches   = D*1000/25.4; % Diameter in inches
rps_max    = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter

%% Electric
tau = 0.2;
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency
mbatt = 7; %[kg]
e0    = 720e3; %[J/kg]
E_batt     = e0*mbatt;  %[J] Total energy of the battery packs

%% Power model

CP0 = 0.034214580122684;
CP1 = -0.002523708791417;
CP2 = 0.116121898742278;
CP3 = -0.248807360672063;



%% Define interpolations
Nlin = 20;  % RPS cases
Tlin = 100; % Steps for time
n_v = linspace(65,rps_max, Nlin); % Interpolate number of cases for RPS
tspan = linspace(0,30,Tlin);      % Interpolate time of integration
y0 = [0 0];                       % Initial conditions. H(0) = 0; V(0) = 0
TOL = 1e-10;                      % Define tolerance

%% Definition of options
optionsolver = odeset('AbsTol',TOL,'RelTol',TOL,'Events',@eventsH);

%% Preallocate solution vectors for speed
H_mat = zeros(Tlin,Nlin);
V_mat = zeros(Tlin,Nlin);
t_mat = zeros(Tlin,Nlin);

P_mat = zeros(Tlin,Nlin);
E_v   = zeros(1,Nlin);
E_dist = zeros(Tlin,Nlin);
E_dist(1,:) = E_batt;

%% Solve system for each rps
for ii = 1:Nlin
    
    [t,y] = ode45(@(t,y)fun(t,y,n_v(ii)),tspan,y0,optionsolver);
    
    % y
    % length(y(:,1))
    
    %% Save each solution colum into a defined matrix
    
    
    H_mat(:,ii) = interp1(linspace(0,1,numel(y(:,1))),y(:,1),linspace(0,1,Tlin)); % For height
    V_mat(:,ii) = interp1(linspace(0,1,numel(y(:,2))),y(:,2),linspace(0,1,Tlin)); % For speed
    t_mat(:,ii) = interp1(linspace(0,1,numel(t)),t,linspace(0,1,Tlin));           % For time
    
    %% Energy Calculations
    
    
    %% Define CP as a vector and not as a function to increase calculation speed within the loop
    P_mat (:,ii) =  (CP3*V_mat(:,ii).^3./(n_v(ii).^3*D^3)+CP2*V_mat(:,ii).^2./(n_v(ii).^2*D^2)+CP1*V_mat(:,ii)./(n_v(ii).*D)+CP0).*rho*N_eng*n_v(ii).^3*D^5;
    
    E_v   (ii)   = trapz(t_mat(:,ii),P_mat (:,ii))/((1-tau)*eta_m);
    
    %% To calculate Energy distribution over time
    for jj = 2:Tlin
        E_dist(jj,ii) = E_batt - trapz(t_mat(1:jj,ii),P_mat(1:jj,ii))/((1-tau)*eta_m);
        
    end
    
end




%% Minimum calculated energy
[minEv,indexminE] = min(E_v);

%% Calculate optimum rps solution

% RPS that give minimum Energy consumption
n_min = n_v(indexminE);


%% Define matrices for plot
nmat = repmat(n_v,Tlin,1);

vect_cc_n = linspace(min(min(nmat)),max(max(nmat)),Nlin)';



%% Plot dynamics

figure(1)
[n_c,h_nmat_c] = contourf([H_mat(:,1:indexminE-1) H_mat(:,indexminE+1:Nlin)],...
    [V_mat(:,1:indexminE-1) V_mat(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'b','LineWidth',1.5);
clabel(n_c,h_nmat_c)
hold on
plot (H_mat(:,indexminE),V_mat(:,indexminE),'r','LineWidth',1.5);

colormap(flipud(colormap('white')))
grid on
Title1 = strcat(['Vertical speed [m/s] vs Height [m] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
title(Title1)
xlabel('Climb height [m] ')
ylabel('V [m/s]')
caxis([min(min(nmat)) max(max(nmat))]); %%%%%%%%%%%%%%%
%legend('E (V,n)','T = D + W constrain','Max RPS','Optimal solution Analytical','Optimal solution Fmincon','Location','northwest')

figure(2)

plot(n_v,E_v,'b','LineWidth',1.8)
grid on
hold on
plot(n_min,minEv,'or','LineWidth',1.2)
xlabel('Engine revolutions [rps] ')
ylabel('Energy consumption [J]')
Title1 = strcat('Energy consumption vs Engine rev. [rps]');
title(Title1)
legend('E(n)',['Opt Energy consumption = ' num2str(minEv/1000) ' kJ.'],'Location','best')

figure(3)

[n_c2,t_nmat_c] = contourf([t_mat(:,1:indexminE-1) t_mat(:,indexminE+1:Nlin)],...
    [V_mat(:,1:indexminE-1) V_mat(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'k','LineWidth',1.5);
clabel(n_c2,t_nmat_c)
colormap(flipud(colormap('white')))
hold on
plot (t_mat(:,indexminE),V_mat(:,indexminE),'r','LineWidth',1.5);
grid on
Title3 = strcat(['Vertical speed [m/s] vs time [s] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
title(Title3)
xlabel('Climb time [s] ')
ylabel('V [m/s]')
caxis([min(min(nmat)) max(max(nmat))]); %%%%%%%%%%%%%%%


figure(4)

[n_c3,ht_nmat_c] = contourf([t_mat(:,1:indexminE-1) t_mat(:,indexminE+1:Nlin)],...
    [H_mat(:,1:indexminE-1) H_mat(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'k','LineWidth',1.5);
clabel(n_c3,ht_nmat_c)
colormap(flipud(colormap('gray')))
hold on
plot (t_mat(:,indexminE),H_mat(:,indexminE),'r','LineWidth',1.5);
grid on
Title3 = strcat(['Climb height [m] vs time [s] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
title(Title3)
xlabel('Climb time [s] ')
ylabel('Climb height [m]')
sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
caxis([min(min(nmat)) 1.3*max(max(nmat))]); %%%%%%%%%%%%%%%


figure(5)
[n_c4,E_nmat_c] = contourf([t_mat(:,1:indexminE-1) t_mat(:,indexminE+1:Nlin)],...
    [E_dist(:,1:indexminE-1) E_dist(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'k','LineWidth',1.5);
clabel(n_c4,E_nmat_c)
colormap(flipud(colormap('white')))
hold on
plot (t_mat(:,indexminE),E_dist(:,indexminE),'r','LineWidth',1.5);
grid on
Title3 = strcat(['Energy storage vs time [s] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
title(Title3)
xlabel('Climb time [s] ')
ylabel('E [J]')
caxis([min(min(nmat)) max(max(nmat))]); %%%%%%%%%%%%%%%
%% Definition of function to be solved

function dy = fun(t,y,n)
% Some parameters are not required within this function, but they are
% included anyway. Will delete in a further version of the code.

global D
global N_eng
global rho

%% Aerodynamic model
alpha_V = -90*pi/180; % [º] Because of Vertical TO


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


% Full model w.r.t angle of attack

CD = @(alpha) p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + p_CD_ac(4)*alpha.^3 ...
    + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);

%% Propulsive Model

CT0 =  0.089050757500461;
CT1 = -0.026659496766654;
CT2 = -0.162620845549081;




%% Loads and weight
mTOW = 21.39; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
%W     = mTOW*g;

S = 0.430279179101; % Reference surface [m^2]




%% Thrust
CT = @(V,n) CT2*V^2/(n^2*D^2)+CT1*V/(n*D)+CT0;

T = @(V,n) CT(V,n)*N_eng*rho*n^2*D^4;




%% Define system of equations

dy = zeros(2,1);
dy(1) = y(2); % dh/dt = V
dy(2) = T(y(2),n)/mTOW - g - 1/2*rho*y(2).^2*S*CD(alpha_V)/mTOW; % m*dV/dt = T-W-D

end


function [position,isterminal,direction] = eventsH(t,y)
global deltaH
position = y(1)-deltaH; % The value that we want to be zero
isterminal = 1;           % Halt integration
direction = 0;            % The zero can be approached from either direction
end


