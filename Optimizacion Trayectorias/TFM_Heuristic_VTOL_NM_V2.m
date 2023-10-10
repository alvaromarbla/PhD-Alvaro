clear all
close all
clc


%h_int = 65;
deltaHtot = 50; %[m]
deltaH_0 = deltaHtot/2;
flagplot = 0;
h_sol  = fzero (@(h_int) Heuristic_VTOL(h_int,flagplot)-deltaHtot,deltaH_0); %Heuristic_VTOL (h_int)

flagplot = 1;

Heuristic_VTOL(h_sol,flagplot);
function H_end = Heuristic_VTOL (h_int,flagplot)

%% Code for solving Min Enegy consumption for a altitude difference for Vertical Take-off

%% Set deltaH (altitude difference)
deltaH1 = h_int;
%% Main parameters definition
D     = 0.8128; % Propeller Diameter [m]
N_eng = 2; % Number of engines
rho   = 1.2133; % Air density [kg/m^3]

%% Define maximum rps


RPMMAX_APC = 150000; % Max RPM of AXI motors
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches*60); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;
%% Electric
tau = 0.2;
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency
mbatt = 3; %[kg]
e0    = 720e3; %[J/kg]
E_batt     = e0*mbatt;  %[J] Total energy of the battery packs

%% Power model


CP0 = 0.0261518307541734;
CP1 = 0.0473735972985378;
CP2 = -0.16267474946046;
CP3 = 0.0247028469343899;

CP = @(V,n) CP3*V.^3./(n^3*D^3)+CP2*V.^2./(n^2*D^2)+CP1*V./(n*D)+CP0;

P = @(V,n) CP(V,n)*N_eng*rho*n^3*D^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model definition done. Set up ode solver.

%% Define interpolations
%Nlin = 20;  % RPS cases
Tlin = 500; % Steps for time
n_v = rps_max*0.75%linspace(50,rps_max, Nlin); % Interpolate number of cases for RPS
tspan = linspace(0,30,Tlin);      % Interpolate time of integration
y0 = [0 0];                       % Initial conditions. H(0) = 0; V(0) = 0
TOL = 1e-10;                      % Define tolerance

%% Definition of options
optionsolver1 = odeset('AbsTol',TOL,'RelTol',TOL,'Events',@eventsH1);



%% Preallocate solution vectors for speed
% H_mat = zeros(Tlin,1);%zeros(Tlin,Nlin);
% V_mat = zeros(Tlin,1);%zeros(Tlin,Nlin);
% t_mat = zeros(Tlin,1);%zeros(Tlin,Nlin);
%
% P_mat = zeros(Tlin,1);%zeros(Tlin,Nlin);
% E_dist = zeros(Tlin,1);%zeros(Tlin,Nlin);
% E_perc = zeros(Tlin,1);%zeros(Tlin,Nlin);

E_dist(1,:) = E_batt;
E_perc(1,:) = 100;

[t1,y1] = ode45(@(t,y)fun(t,y,n_v,D,N_eng,rho),tspan,y0,optionsolver1);

% Create vector for n_vec for first half
n_vec_1 =  n_v*ones(1,length(t1));
%% First segment is solved. Now solve for the second!

% Extract initial conditions from previous segment.


y0_2(1) = y1(end,1);
y0_2(2) = y1(end,2);
V_break = y1(end,2);
% Define new revolutions

n_Hover = 45.7558;  % Hovering
n_v = n_Hover*0.9;
optionsolver2 = odeset('AbsTol',TOL,'RelTol',TOL,'Events',@eventsH2);

[t2,y2] = ode45(@(t,y)fun(t,y,n_v,D,N_eng,rho),tspan,y0_2,optionsolver2);
t2 = t1(end) + t2;

t_sol = [t1; t2(2:end)];
y_sol = [y1 ; y2(2:end,:)];
%% Save each solution colum into a defined matrix


H_vec = interp1(linspace(0,1,numel(y_sol(:,1))),y_sol(:,1),linspace(0,1,Tlin)); % For height
V_vec = interp1(linspace(0,1,numel(y_sol(:,2))),y_sol(:,2),linspace(0,1,Tlin)); % For speed
t_vec = interp1(linspace(0,1,numel(t_sol)),t_sol,linspace(0,1,Tlin));           % For time
n_vec_2 =  n_v*ones(1,length(t2)-length(t1));

%% Energy Calculations


P_break = CP(V_break,n_v).*rho*N_eng*n_v.^3*D^5;

P_vec  = CP(V_vec,n_v).*rho*N_eng*n_v.^3*D^5;



%E_v      = trapz(t_vec,P_vec)/((1-tau)*eta_m);

%% To calculate Energy distribution over time
for jj = 2:Tlin
    E_dist(jj) = E_batt - trapz(t_vec(1:jj),P_vec(1:jj))/(eta_m);
    E_perc(jj) = E_dist(jj)/E_batt *100;

    if jj == length(t1) 
        E_break = E_dist(jj);
        E_perc_break = E_perc(jj);
    end
end
E_dist(length(t1));
E_break;


H_end = H_vec(end)

if flagplot ==1

figure(1)
plot(t_vec,H_vec,'LineWidth',1.5)
grid on
Title1 = strcat(['Time [s] vs Height [m]. Hover break height at ' num2str(h_int) ' m.']);
title(Title1)
hold on 
plot (t1(end),h_int,'or','LineWidth',1.5)
xlabel('t [s]')
ylabel('h [m]')
legend('Height profile','Break point','Location','best')


figure(2)
plot(t_vec,V_vec,'LineWidth',1.5)
grid on
Title1 = strcat(['Time [s] vs Speed [m]. Hover break speed at ' num2str(V_break) ' m/s.']);
title(Title1)
hold on 
plot (t1(end),V_break,'or','LineWidth',1.5)
xlabel('t [s]')
ylabel('V [m/s]')
legend('Speed profile','Break point','Location','northeast')


figure(3)
plot(t_vec,P_vec,'LineWidth',1.5)
grid on
Title1 = strcat(['Time [s] vs Power [W]. Hover break Power at ' num2str(P_break) ' W.']);
title(Title1)
hold on 
plot (t1(end),P_break,'or','LineWidth',1.5)
xlabel('t [s]')
ylabel('P [W]')
legend('Power profile','Break point','Location','northeast')

figure(4)
plot(t_vec,E_dist,'LineWidth',1.5)
grid on
Title1 = strcat(['Time [s] vs Energy [J]. Hover break Energy at ' num2str(E_break*10e-3) ' kJ.']);
title(Title1)
hold on 
plot (t1(end),E_dist(length(t1))*0.997,'or','LineWidth',1.5)%*0.99935
xlabel('t [s]')
ylabel('E [J]')
legend('Energy profile','Break point','Location','northeast')

figure(5)
plot(t_vec,E_perc,'LineWidth',1.5)
grid on
Title1 = strcat(['Time [s] vs Energy discharge perc. [%]. Hover break Energy at ' num2str(E_perc_break) ' %.']);
title(Title1)
hold on 
plot (t1(end),E_perc(length(t1))*0.997,'or','LineWidth',1.5)
xlabel('t [s]')
ylabel('E [%]')
legend('Energy percentage profile','Break point','Location','northeast')


figure(6)
plot(H_vec,V_vec,'LineWidth',1.5)
grid on
Title1 = strcat(['Height [m] vs Speed [m/s]. Hover break Speed at ' num2str(V_break) ' m/s.']);
title(Title1)
hold on 
plot (h_int,V_break,'or','LineWidth',1.5)
xlabel('h [m]')
ylabel('V [m/s]')
legend('Speed profile','Break point','Location','northeast')

figure(7)
plot(H_vec,P_vec,'LineWidth',1.5)
grid on
Title1 = strcat(['Height [m] vs Power [W]. Hover break Power at ' num2str(P_break) ' W.']);
title(Title1)
hold on 
plot (h_int,P_break,'or','LineWidth',1.5)
xlabel('h [m]')
ylabel('P [W]')
legend('Power profile','Break point','Location','best')

n_vec = [n_vec_1, n_vec_2];
n_vec = interp1(linspace(0,1,numel(n_vec)),n_vec,linspace(0,1,Tlin));
Iterant_mission = [t_vec; H_vec;V_vec;E_dist; n_vec]';

save Iterant_mission.mat Iterant_mission

end
% %% Minimum calculated energy
% [minEv,indexminE] = min(E_v);

%% Calculate optimum rps solution

% % RPS that give minimum Energy consumption
% n_min = n_v(indexminE);


% %% Define matrices for plot
% nmat = repmat(n_v,Tlin,1);
%
% vect_cc_n = linspace(min(min(nmat)),max(max(nmat)),Nlin)';



% %% Plot dynamics
%
% figure(1)
% [n_c,h_nmat_c] = contourf([H_mat(:,1:indexminE-1) H_mat(:,indexminE+1:Nlin)],...
%     [V_mat(:,1:indexminE-1) V_mat(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'b','LineWidth',1.5);
% clabel(n_c,h_nmat_c)
% hold on
% plot (H_mat(:,indexminE),V_mat(:,indexminE),'r','LineWidth',1.5);
%
% colormap(flipud(colormap('white')))
% grid on
% Title1 = strcat(['Vertical speed [m/s] vs Height [m] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
% title(Title1)
% xlabel('Climb height [m] ')
% ylabel('V [m/s]')
% caxis([min(min(nmat)) max(max(nmat))]); %%%%%%%%%%%%%%%
% %legend('E (V,n)','T = D + W constrain','Max RPS','Optimal solution Analytical','Optimal solution Fmincon','Location','northwest')
%
% figure(2)
%
% plot(n_v,E_v,'b','LineWidth',1.8)
% grid on
% hold on
% xline(rps_max,'--b','LineWidth',1.5)
% plot(n_min,minEv,'or','LineWidth',1.2)
% xlabel('Engine revolutions [rps] ')
% ylabel('Energy consumption [J]')
% Title1 = strcat('Energy consumption vs Engine rev. [rps]');
% title(Title1)
% legend('E(n)',['Opt Energy consumption = ' num2str(minEv/1000) ' kJ.'],'Location','best')
%
% figure(3)
%
% [n_c2,t_nmat_c] = contourf([t_mat(:,1:indexminE-1) t_mat(:,indexminE+1:Nlin)],...
%     [V_mat(:,1:indexminE-1) V_mat(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'k','LineWidth',1.5);
% clabel(n_c2,t_nmat_c)
% colormap(flipud(colormap('white')))
% hold on
% plot (t_mat(:,indexminE),V_mat(:,indexminE),'r','LineWidth',1.5);
% grid on
% Title3 = strcat(['Vertical speed [m/s] vs time [s] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
% title(Title3)
% xlabel('Climb time [s] ')
% ylabel('V [m/s]')
% caxis([min(min(nmat)) max(max(nmat))]); %%%%%%%%%%%%%%%
%
%
% figure(4)
%
% [n_c3,ht_nmat_c] = contourf([t_mat(:,1:indexminE-1) t_mat(:,indexminE+1:Nlin)],...
%     [H_mat(:,1:indexminE-1) H_mat(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'k','LineWidth',1.5);
% clabel(n_c3,ht_nmat_c)
% colormap(flipud(colormap('gray')))
% hold on
% plot (t_mat(:,indexminE),H_mat(:,indexminE),'r','LineWidth',1.5);
% grid on
% Title3 = strcat(['Climb height [m] vs time [s] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
% title(Title3)
% xlabel('Climb time [s] ')
% ylabel('Climb height [m]')
% sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
% sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
% caxis([min(min(nmat)) 1.3*max(max(nmat))]); %%%%%%%%%%%%%%%
%
%
% figure(5)
% [n_c4,E_nmat_c] = contourf([t_mat(:,1:indexminE-1) t_mat(:,indexminE+1:Nlin)],...
%     [E_dist(:,1:indexminE-1) E_dist(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'k','LineWidth',1.5);
% clabel(n_c4,E_nmat_c)
% colormap(flipud(colormap('white')))
% hold on
% plot (t_mat(:,indexminE),E_dist(:,indexminE),'r','LineWidth',1.5);
% grid on
% Title3 = strcat(['Energy storage vs time [s] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
% title(Title3)
% xlabel('Climb time [s] ')
% ylabel('E [J]')
% caxis([min(min(nmat)) max(max(nmat))]); %%%%%%%%%%%%%%%
%
% figure(6)
% [n_c4,Eper_nmat_c] = contourf([t_mat(:,1:indexminE-1) t_mat(:,indexminE+1:Nlin)],...
%     [E_perc(:,1:indexminE-1) E_perc(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'k','LineWidth',1.5);
% clabel(n_c4,Eper_nmat_c)
% colormap(flipud(colormap('white')))
% hold on
% plot (t_mat(:,indexminE),E_perc(:,indexminE),'r','LineWidth',1.5);
% grid on
% Title3 = strcat(['Energy percentage [%] vs time [s] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
% title(Title3)
% xlabel('Climb time [s] ')
% ylabel('E [%]')
% caxis([min(min(nmat)) max(max(nmat))]); %%%%%%%%%%%%%%%
%
%
% figure(7)
%
% [n_c2,P_nmat_c] = contourf([H_mat(:,1:indexminE-1) H_mat(:,indexminE+1:Nlin)],...
%     [P_mat(:,1:indexminE-1) P_mat(:,indexminE+1:Nlin)],[nmat(:,1:indexminE-1) nmat(:,indexminE+1:Nlin)],[vect_cc_n(1:indexminE-1); vect_cc_n(indexminE+1:Nlin)],'k','LineWidth',1.5);
% clabel(n_c2,P_nmat_c)
% colormap(flipud(colormap('white')))
% hold on
% plot (H_mat(:,indexminE),P_mat(:,indexminE),'r','LineWidth',1.5);
% %yline(Pmax,'-.b','LineWidth',2);
% grid on
% Title3 = strcat(['Power consumption [W] vs Height [m] & Engine rev. [rps]. Optimum at ' num2str(n_min) ' rps.']);
% title(Title3)
% xlabel('Climb height [m] ')
% ylabel('P [W]')
% caxis([min(min(nmat)) max(max(nmat))]); %%%%%%%%%%%%%%%


%% Definition of function to be solved

    function dy = fun(t,y,n,D,N_eng,rho)
        % Some parameters are not required within this function, but they are
        % included anyway. Will delete in a further version of the code.


        %% Aerodynamic model
        alpha_V = -90*pi/180; % [º] Because of Vertical TO


        p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];


        % Full model w.r.t angle of attack

        CD = @(alpha) p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + p_CD_ac(4)*alpha.^3 ...
            + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);





        %% Loads and weight
        mTOW = 16.6; % Maximum T-O Mass [kg]
        g    = 9.81; % Gravity [m/s^2]
        %W     = mTOW*g;

        S = 0.5004; % Reference surface [m^2]

        %% Thrust
        CT0 =  0.0735880531010883;
        CT1 = -0.0311758018412727;
        CT2 = -0.249744726429543;

        CT = @(V,n) CT2*V^2/(n^2*D^2)+CT1*V/(n*D)+CT0;

        T = @(V,n) CT(V,n)*N_eng*rho*n^2*D^4;


        %% Define system of equations

        dy = zeros(2,1);
        dy(1) = y(2); % dh/dt = V
        dy(2) = T(y(2),n)/mTOW - g - 1/2*rho*y(2).^2*S*CD(alpha_V)/mTOW; % m*dV/dt = T-W-D

    end


    function [position,isterminal,direction] = eventsH1(t,y)
        position = y(1)-deltaH1; % The value that we want to be zero
        isterminal = 1;           % Halt integration
        direction = 1;            % The zero can be approached from either direction
    end

function [position,isterminal,direction] = eventsH2(t,y)
        position = y(2); % The value that we want to be zero
        isterminal = 1;           % Halt integration
        direction = -1;            % The zero can be approached from either direction
    end

end

