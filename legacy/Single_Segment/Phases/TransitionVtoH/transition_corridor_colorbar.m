%function transition_corridor_colorbar()
% This script calculates the transition corridor.

colormap_flag = 1;
% Lower limit is calculated from condition of max operative angle of
% attack, solving n.

% Higher limit is calculated from max n and solving for alpha. 
params = Transition_Config();
eps_min = -5*pi/180;
eps_max = pi/2;
epsilon_vec = linspace(eps_min,eps_max,200);


%% Lower Limit

alpha = params.alpha_max_ope;
V_sol_low = zeros(1, length(epsilon_vec));
n_sol_low = zeros(1, length(epsilon_vec));
exitflag_low = zeros(1, length(epsilon_vec));
fval_low    = zeros(1, length(epsilon_vec));

for ii = 1:length(epsilon_vec)

    % Reuse models from descent
    T = @(X) CT_Model(X(1), X(2), alpha, epsilon_vec(ii), params);
    L = @(X) CL_Model(alpha) * 0.5*params.rho*X(1)^2*params.S_ref;
    D = @(X) CD_Model(alpha) * 0.5*params.rho*X(1)^2*params.S_ref;
    fun = @(X)[(T(X)*cos(epsilon_vec(ii)+alpha) - D(X))/(params.mass*params.g);   % Axial force balance
        (L(X) + T(X)*sin(epsilon_vec(ii)+alpha))/(params.mass*params.g) - 1];  % Normal force balance

    [X_sol,fval, exitflag_low(ii)] = fsolve( fun, [19,20]);
    fval_low (ii)= norm(fval);
    V_sol_low(ii) = X_sol(1);
    n_sol_low(ii)  = X_sol(2);
end

T_debug = CT_Model(V_sol_low, n_sol_low, -5*pi/180, epsilon_vec(ii), params);
D_debug = CD_Model(alpha)*0.5*params.rho*V_sol_low.^2*params.S_ref;
L_debug = CL_Model(alpha)*0.5*params.rho*V_sol_low.^2*params.S_ref;
J_debug = V_sol_low./(params.prop.diameter*n_sol_low);
matriz_Res = [180/pi*epsilon_vec' V_sol_low', n_sol_low', T_debug' D_debug' J_debug' L_debug' exitflag_low' fval_low']

disp('End of Low')
%% Higher Limit

alpha = -5*pi/180;
V_sol_hi = zeros(1, length(epsilon_vec));
n_sol_hi = zeros(1, length(epsilon_vec));

for ii = 1:length(epsilon_vec)

    % Reuse models from descent
    T = @(X) CT_Model(X(1), X(2), alpha, epsilon_vec(ii), params);
    L = @(X) CL_Model(alpha) * 0.5*params.rho*X(1)^2*params.S_ref;
    D = @(X) CD_Model(alpha) * 0.5*params.rho*X(1)^2*params.S_ref;
    fun = @(X)[T(X)*cos(epsilon_vec(ii)+alpha) - D(X);   % Axial force balance
        L(X) + T(X)*sin(epsilon_vec(ii)+alpha) - params.mass*params.g];  % Normal force balance

    X_sol = fsolve( fun, [30,5*pi/180]);
    V_sol_hi(ii) = X_sol(1);
    n_sol_hi(ii) = X_sol(2);
end

%% Refactor curves

kk = find (V_sol_low<1e-5);
first_zero = kk(1);

low_eps = linspace(V_sol_low(1),V_sol_hi(1),50);
high_eps = linspace(0,V_sol_hi(end),50);

vert_vel = linspace(V_sol_low(first_zero),eps_max);
%% Plot

hFIG = figure;


plot( V_sol_low(1:first_zero),epsilon_vec(1:first_zero)*180/pi , 'LineWidth', 1.5, 'Color', 'r')
xlim ([0 50])
ylim ([-5 90])
hold on
plot( V_sol_hi,epsilon_vec*180/pi, 'LineWidth', 1.5 ,'Color', 'r')
plot(high_eps,90*ones(1,length(high_eps)),'LineWidth', 1.5 ,'Color', 'k')
plot(zeros(1,length(vert_vel)),linspace(epsilon_vec(first_zero)*180/pi,90,length(vert_vel)),'LineWidth', 1.5 ,'Color', 'k')
plot(low_eps,eps_min*ones(1,length(low_eps))*180/pi,'LineWidth', 1.5 ,'Color', 'k')

if colormap_flag ==1
    scatter(V_sol_hi,epsilon_vec*180/pi,20,n_sol_hi,'filled')
    scatter(V_sol_low,epsilon_vec*180/pi,20,n_sol_low,'filled')

    colormap('jet')
    colorbar
end

X_fill = [V_sol_low(1:first_zero),...
    zeros(1,length(vert_vel)),...
    high_eps,...
    fliplr(V_sol_hi),...
    fliplr(low_eps)];

Y_fill = [epsilon_vec(1:first_zero)*180/pi,...
    linspace(epsilon_vec(first_zero)*180/pi,90,length(vert_vel)),...
    90*ones(1,length(high_eps)),...
    fliplr(epsilon_vec*180/pi)...
    fliplr(180/pi*eps_min*ones(1,length(low_eps)))];

% Create ylabel
ylabel('$\varepsilon$ [deg]');

% Create xlabel
xlabel('V [m/s]');

grid on
% Create title
title('Transition corridor');


fill(X_fill,Y_fill,[0.7,0.7,0.7],'FaceAlpha',.2)

fname = "transition_corridor";


picturewidth = 20;
hw_ratio = 0.65;
set(findall(hFIG,'-property','FontSize'),'FontSize',16)
text(5,30,"Low Speed Limit",'FontSize',12);
text(35,75,"High Speed Limit",'FontSize',12);
set(findall(hFIG,'-property','Box'),'Box','off') % optional
set(findall(hFIG,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hFIG,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hFIG,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hFIG,'Position');
set(hFIG,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hFIG,fname,'-dpdf','-vector','-fillpage')

