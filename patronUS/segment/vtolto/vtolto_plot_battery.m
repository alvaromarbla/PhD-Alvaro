function vtolto_plot_battery(res, params, name)
% Battery page of a trajectory segment (VTOL take-off, time-dependent
% cruise) with the Shepherd model (res from
% unpack_ocp_solution): bus power, current, terminal voltage against the
% cut-off and open-circuit voltage, state of charge against the
% ideal-battery energy fraction, and voltage against charge drawn, both
% along the take-off and over a full discharge at constant power (the
% classic discharge curve, with the take-off's operating points on it).
% name: output PDF in Images/ (default 'vtolto_battery')

if nargin < 3, name = 'vtolto_battery'; end
if ~isfield(res, 'SOC'), return; end

t = res.t;
b = params.batt;
U_oc = b.E0 + b.A*exp(-b.B*(1 - res.SOC));

battFig = figure('Name', ['Battery: ' name], 'Position', [100 100 1500 700]);
if t(end) > 120                                % Long segments (cruise) in minutes
    t = t/60; t_label = '$t$ [min]';
else
    t_label = '$t$ [s]';
end
drawn = 100*(1 - res.SOC);                     % Charge drawn [% of capacity]

subplot(2,3,1); hold on; grid on;
plot(t, res.dbg.P_bus/1e3, 'LineWidth', 1.5);
P_max = U_oc.^2./(4*(b.R + b.K./res.SOC));
plot(t, P_max/1e3, 'r--', 'LineWidth', 1.2);
xlabel(t_label); ylabel('[kW]'); title('Bus power'); legend('$P_{bus}$', 'max deliverable', 'Location', 'best');

subplot(2,3,2); hold on; grid on;
plot(t, res.dbg.I_batt, 'LineWidth', 1.5);
xlabel(t_label); ylabel('[A]'); title('Battery current');

subplot(2,3,3); hold on; grid on;
plot(t, res.dbg.U_batt, 'LineWidth', 1.5);
plot(t, U_oc, ':', 'LineWidth', 1.2);
yline(b.V_cut, 'r--');
xlabel(t_label); ylabel('[V]'); title('Terminal voltage'); legend('$U$', '$U_{oc}$', 'cut-off', 'Location', 'best');

subplot(2,3,4); hold on; grid on;
plot(t, drawn, 'LineWidth', 1.5);
plot(t, 100*res.E/params.E_batt, '--', 'LineWidth', 1.5);
xlabel(t_label); ylabel('[\%]'); title('Charge and energy drawn');
legend('$1 - SOC$ (Shepherd)', '$E/E_{batt}$ (ideal)', 'Location', 'best');

% Voltage against charge drawn along the take-off: the sag follows the current
subplot(2,3,5); hold on; grid on;
plot(drawn, res.dbg.U_batt, '-o', 'LineWidth', 1.5, 'MarkerSize', 3);
plot(drawn, U_oc, ':', 'LineWidth', 1.2);
yline(b.V_cut, 'r--');
xlabel('charge drawn [\%]'); ylabel('[V]'); title('Voltage vs charge drawn');
legend('$U$', '$U_{oc}$', 'cut-off', 'Location', 'best');

% Full discharge at constant power: U(SOC) from the same Shepherd law, at
% the take-off's final hover power and at its peak power
subplot(2,3,6); hold on; grid on;
SOC_grid = linspace(1, 0.02, 400);
U_oc_grid = b.E0 + b.A*exp(-b.B*(1 - SOC_grid));
plot(100*(1 - SOC_grid), U_oc_grid, ':', 'LineWidth', 1.2, 'DisplayName', '$U_{oc}$');
P_levels = [res.dbg.P_bus(end), max(res.dbg.P_bus)];
P_names  = {'$t_f$', 'peak'};
for i = 1:2
    op = battery_shepherd(P_levels(i), SOC_grid, b);
    U  = op.U;
    U(op.disc < 0) = NaN;                       % Power not deliverable
    plot(100*(1 - SOC_grid), U, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('$U$ at %s, %.1f kW', P_names{i}, P_levels(i)/1e3));
end
plot(drawn, res.dbg.U_batt, 'k.', 'MarkerSize', 10, 'DisplayName', 'trajectory');
yline(b.V_cut, 'r--', 'HandleVisibility', 'off');
xlabel('charge drawn [\%]'); ylabel('[V]'); title('Full discharge at constant power');
legend('Location', 'southwest');

picturewidth = 60;
hw_ratio = 0.5;
set(findall(battFig,'-property','FontSize'),'FontSize',16)
set(findall(battFig,'-property','Interpreter'),'Interpreter','latex')
set(findall(battFig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(battFig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(battFig,'Position');
set(battFig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
targetFolder = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'Images');
if ~isfolder(targetFolder), mkdir(targetFolder); end
print(battFig, fullfile(targetFolder, name),'-dpdf','-vector','-fillpage')

end
