function vtolto_plot_results(res, params)
% Results page of the VTOL take-off (res from unpack_ocp_solution)

t = res.t;

vtoltoFig = figure('Name', 'VTOL Take-Off Results', 'Position', [100 100 1100 800]);

subplot(3,2,1); hold on; grid on;
plot(res.x, res.h, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
axis equal; xlabel('$x$ [m]'); ylabel('$h$ [m]'); title('Trajectory');

subplot(3,2,2); hold on; grid on;
plot(t, res.vx, 'LineWidth', 1.5); plot(t, res.vh, 'LineWidth', 1.5); plot(t, res.dbg.V, 'k--', 'LineWidth', 1.5);
xlabel('$t$ [s]'); ylabel('[m/s]'); title('Velocities'); legend('$v_x$', '$v_h$', '$V$', 'Location', 'best');

subplot(3,2,3); hold on; grid on;
plot(t, rad2deg(res.theta), 'LineWidth', 1.5); plot(t, rad2deg(res.dbg.alpha), 'LineWidth', 1.5); plot(t, rad2deg(res.dbg.gamma), 'LineWidth', 1.5);
xlabel('$t$ [s]'); ylabel('[deg]'); title('Attitude'); legend('$\theta$', '$\alpha$', '$\gamma$', 'Location', 'best');

subplot(3,2,4); hold on; grid on;
plot(t, rad2deg(res.epsilon1), 'LineWidth', 1.5); plot(t, rad2deg(res.epsilon2), 'LineWidth', 1.5);
xlabel('$t$ [s]'); ylabel('[deg]'); title('Tilt angles'); legend('$\epsilon_1$', '$\epsilon_2$', 'Location', 'best');

subplot(3,2,5); hold on; grid on;
plot(t, res.n1, 'LineWidth', 1.5); plot(t, res.n2, 'LineWidth', 1.5);
yline(params.prop.max_rps, 'r--');
xlabel('$t$ [s]'); ylabel('[rps]'); title('Propeller speed'); legend('$n_1$', '$n_2$', 'Location', 'best');

subplot(3,2,6); hold on; grid on;
yyaxis left;  plot(t, (res.dbg.P1 + res.dbg.P2)/1e3, 'LineWidth', 1.5); ylabel('$P$ [kW]');
yyaxis right; plot(t, 100*res.E/params.E_batt, 'LineWidth', 1.5); ylabel('$E/E_{batt}$ [\%]');
xlabel('$t$ [s]'); title('Power and energy');

picturewidth = 50;
hw_ratio = 0.75;
set(findall(vtoltoFig,'-property','FontSize'),'FontSize',16)
set(findall(vtoltoFig,'-property','Interpreter'),'Interpreter','latex')
set(findall(vtoltoFig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(vtoltoFig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(vtoltoFig,'Position');
set(vtoltoFig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
targetFolder = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'Images');
if ~isfolder(targetFolder), mkdir(targetFolder); end
print(vtoltoFig, fullfile(targetFolder, 'vtolto_results'),'-dpdf','-vector','-fillpage')

end
