function plot_trajectory_results(fh, params, name)
% Results page of a trajectory segment over time (fh from flight_history,
% so the same page serves the VTOL take-off, the time-dependent cruise and,
% with concatenated histories, a multiphase trajectory): path, velocities,
% attitude, tilt angles, rotor speeds, power and energy.
% name: output PDF in Images/ (default [fh.segment '_results'])

if nargin < 3, name = [fh.segment '_results']; end

[t, t_label] = time_axis(fh.t);

fig = figure('Name', ['Results: ' name], 'Position', [100 100 1100 800]);

% Path: x-h, or range over time when the altitude does not change (cruise)
subplot(3,2,1); hold on; grid on;
long = max(fh.x) - min(fh.x) > 2e3;
if range(fh.h) < 1
    plot(t, fh.x/1e3, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
    xlabel(t_label); ylabel('$x$ [km]');
    title(sprintf('Range ($h$ = %.0f m)', fh.h(1)));
elseif long
    plot(fh.x/1e3, fh.h, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
    xlabel('$x$ [km]'); ylabel('$h$ [m]'); title('Trajectory');
else
    plot(fh.x, fh.h, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
    axis equal; xlabel('$x$ [m]'); ylabel('$h$ [m]'); title('Trajectory');
end

subplot(3,2,2); hold on; grid on;
plot(t, fh.vx, 'LineWidth', 1.5); plot(t, fh.vh, 'LineWidth', 1.5); plot(t, fh.V, 'k--', 'LineWidth', 1.5);
xlabel(t_label); ylabel('[m/s]'); title('Velocities'); legend('$v_x$', '$v_h$', '$V$', 'Location', 'best');

subplot(3,2,3); hold on; grid on;
plot(t, rad2deg(fh.theta), 'LineWidth', 1.5); plot(t, rad2deg(fh.alpha), 'LineWidth', 1.5); plot(t, rad2deg(fh.gamma), 'LineWidth', 1.5);
names = {'$\theta$', '$\alpha$', '$\gamma$'};
if ~all(isnan(fh.deltae))
    plot(t, rad2deg(fh.deltae), 'LineWidth', 1.5); names{end+1} = '$\delta_e$';
end
xlabel(t_label); ylabel('[deg]'); title('Attitude'); legend(names, 'Location', 'best');

subplot(3,2,4); hold on; grid on;
plot(t, rad2deg(fh.epsilon1), 'LineWidth', 1.5); plot(t, rad2deg(fh.epsilon2), 'LineWidth', 1.5);
xlabel(t_label); ylabel('[deg]'); title('Tilt angles'); legend('$\epsilon_1$', '$\epsilon_2$', 'Location', 'best');

subplot(3,2,5); hold on; grid on;
plot(t, fh.n1, 'LineWidth', 1.5); plot(t, fh.n2, 'LineWidth', 1.5);
yline(params.prop.max_rps, 'r--');
xlabel(t_label); ylabel('[rps]'); title('Propeller speed'); legend('$n_1$', '$n_2$', 'Location', 'best');

subplot(3,2,6); hold on; grid on;
yyaxis left;  plot(t, fh.P_bus/1e3, 'LineWidth', 1.5); ylabel('$P_{bus}$ [kW]');
yyaxis right; plot(t, 100*fh.E/params.E_batt, 'LineWidth', 1.5); ylabel('[\%]');
names = {'$P_{bus}$', '$E/E_{batt}$'};
if ~all(isnan(fh.SOC))
    plot(t, 100*(1 - fh.SOC), '--', 'LineWidth', 1.5); names{end+1} = '$1 - SOC$';
end
xlabel(t_label); title('Power and energy'); legend(names, 'Location', 'best');

picturewidth = 50;
hw_ratio = 0.75;
set(findall(fig,'-property','FontSize'),'FontSize',16)
set(findall(fig,'-property','Interpreter'),'Interpreter','latex')
set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(fig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
targetFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Images');
if ~isfolder(targetFolder), mkdir(targetFolder); end
print(fig, fullfile(targetFolder, name),'-dpdf','-vector','-fillpage')

end

function [t, label] = time_axis(t)
% Seconds, or minutes for segments longer than 2 min
if t(end) > 120
    t = t/60; label = '$t$ [min]';
else
    label = '$t$ [s]';
end
end
