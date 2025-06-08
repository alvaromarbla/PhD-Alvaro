C_val = 20.039;
C_cost_vec = linspace(0,1.5*C_val,50);

for ii = 1:length(C_cost_vec)
    [params, bounds] = Climb_Config(C_cost_vec(ii));
    sol = Main_Climb_Optimizer(C_cost_vec(ii));
    % Extract optimized parameters
    V = sol(1);             % Airspeed [m/s]
    gamma = sol(2);         % Flight path angle [rad]
    alpha = sol(3);         % Angle of attack [rad]
    epsilon = sol(4);       % Thrust vector angle [rad]
    n = params.prop.max_rps/2;             % Motor RPM
    CP = CP_Model(V, n, alpha, epsilon, params);
    t_climb = params.deltaH / (V*sin(gamma));
    xf(ii) =  params.deltaH * abs(cot(gamma));
    Energy(ii) = params.prop.num_engines * CP * params.rho * n^3 * params.prop.diameter^5 * t_climb/params.prop.eff;
end

hFIG = figure;
fname = "Pareto_Climb";

diff_E = Energy(2:end)-Energy(1:end-1);
diff_xf = xf(2:end)-xf(1:end-1);
Tol_der = 1;
ii_C = find(min((diff_E./diff_xf-C_val).^2)==(diff_E./diff_xf-C_val).^2);
plot(xf,Energy,'DisplayName','Pareto front','LineWidth',1.5)
%axis([0 max(xf) 0 max(Energy)])
hold on 
sl_diff = diff_E(ii_C)/diff_xf(ii_C);
plot (xf,sl_diff.*(xf-xf(ii_C))+Energy(ii_C),'--k','DisplayName','Tangent $C_{Cruise}$','LineWidth',1.5)
plot(xf(ii_C),Energy(ii_C),'DisplayName','$C_{Cruise}$','Marker','o',...
    'LineWidth',1.5,...
    'LineStyle','none',...
    'Color',[1 0 0]);
line([xf(ii_C) xf(ii_C)],[1.5e5 Energy(ii_C)],'Color','k');
text(xf(ii_C),1.65e5,'C = $C_{Cruise}$','Rotation',0,'VerticalAlignment','bottom','HorizontalAlignment','left');
line([xf(1) xf(1)],[1.3e5 Energy(1)],'Color','k');
text(xf(1),1.65e5,'C = 0','Rotation',0,'VerticalAlignment','bottom','HorizontalAlignment','left');
line([xf(end) xf(end)],[1.3e5 Energy(end)],'Color','k');
text(xf(end),1.60e5,'C = 3/2$C_{Cruise}$','Rotation',0,'VerticalAlignment','bottom','HorizontalAlignment','right');
% Create ylabel
ylabel('E [J]');

% Create xlabel
xlabel('$x_f$ [m]');

% Create title
title('Trade-off study for $x_f$ vs $E$');

legend('Convex Envelope of the Pareto Front','Tangent Line','Location','northwest')


grid on

picturewidth = 20;
hw_ratio = 0.65;
set(findall(hFIG,'-property','FontSize'),'FontSize',16)
set(findall(hFIG,'-property','Box'),'Box','off') % optional
set(findall(hFIG,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hFIG,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hFIG,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hFIG,'Position');
set(hFIG,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hFIG,fname,'-dpdf','-vector','-fillpage')