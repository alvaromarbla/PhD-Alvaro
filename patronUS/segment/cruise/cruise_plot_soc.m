function cruise_plot_soc(res, params, guess, name)
% Results page of the time-dependent cruise with the Shepherd battery (res
% from unpack_ocp_solution): range, airspeed, attitude, tilt angles,
% rotor speeds and bus energy per distance along the cruise, against the
% steady optimum held constant (guess, dashed).
% name: output PDF in Images/ (default 'cruise_soc_results')

if nargin < 4, name = 'cruise_soc_results'; end

t_min  = res.t/60;
tg_min = guess.t/60;
red  = [0.85 0.33 0.10];
blue = [0 0.45 0.74];

fig = figure('Name', 'Cruise with Shepherd battery', 'Position', [100 100 1500 700]);

subplot(2,3,1); hold on; grid on;
plot(t_min, res.x/1e3, 'LineWidth', 1.5);
plot(tg_min, guess.x/1e3, '--', 'LineWidth', 1.2);
xlabel('$t$ [min]'); ylabel('[km]'); title('Range');
legend('optimal', 'steady trim held', 'Location', 'best');

subplot(2,3,2); hold on; grid on;
plot(t_min, res.V, 'LineWidth', 1.5);
plot(tg_min, guess.V, '--', 'LineWidth', 1.2);
xlabel('$t$ [min]'); ylabel('[m/s]'); title('Airspeed');

subplot(2,3,3); hold on; grid on;
plot(t_min, rad2deg(res.theta), 'LineWidth', 1.5, 'DisplayName', '$\theta$');
plot(t_min, rad2deg(res.dbg.alpha), ':', 'LineWidth', 1.5, 'DisplayName', '$\alpha$');
if isfield(res, 'deltae')
    plot(t_min, rad2deg(res.deltae), '-.', 'LineWidth', 1.5, 'DisplayName', '$\delta_e$');
end
xlabel('$t$ [min]'); ylabel('[deg]'); title('Attitude'); legend('Location', 'best');

subplot(2,3,4); hold on; grid on;
plot(t_min, rad2deg(res.epsilon1), 'Color', blue, 'LineWidth', 1.5);
plot(t_min, rad2deg(res.epsilon2), 'Color', red,  'LineWidth', 1.5);
xlabel('$t$ [min]'); ylabel('[deg]'); title('Tilt angles');
legend('$\epsilon_1$', '$\epsilon_2$', 'Location', 'best');

subplot(2,3,5); hold on; grid on;
plot(t_min, res.n1, 'Color', blue, 'LineWidth', 1.5);
plot(t_min, res.n2, 'Color', red,  'LineWidth', 1.5);
xlabel('$t$ [min]'); ylabel('[rps]'); title('Rotor speeds');
legend('$n_1$', '$n_2$', 'Location', 'best');

% Bus energy per distance P_bus/V [Wh/km] and charge per distance I/V [Ah/km]
subplot(2,3,6); hold on; grid on;
yyaxis left
plot(t_min, res.dbg.P_bus./res.V/3.6, 'LineWidth', 1.5);
ylabel('$P_{bus}/V$ [Wh/km]');
yyaxis right
plot(t_min, res.dbg.I_batt./res.V/3.6, 'LineWidth', 1.5);
ylabel('$I/V$ [Ah/km]');
xlabel('$t$ [min]'); title('Energy and charge per distance');

picturewidth = 60;
hw_ratio = 0.5;
set(findall(fig,'-property','FontSize'),'FontSize',16)
set(findall(fig,'-property','Interpreter'),'Interpreter','latex')
set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(fig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
targetFolder = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'Images');
if ~isfolder(targetFolder), mkdir(targetFolder); end
print(fig, fullfile(targetFolder, name),'-dpdf','-vector','-fillpage')

end
