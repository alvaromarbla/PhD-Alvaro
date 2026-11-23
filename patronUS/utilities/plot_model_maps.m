function plot_model_maps(op, models, name, t)
% Operating points of a solution on the aerodynamic and propeller models
% (2Wings): where the trajectory, or the steady point, sits in each map.
%   Top:    CT(J, phi) and CP(J, phi) maps with the path of rotor 1 (o) and
%           rotor 2 (s); the bold contour is CT = 0 / CP = 0.
%   Bottom: CL and CD of the wing (alpha), tail (alpha - xi_dw), fuselage
%           (alpha) and nacelles (alpha_nac1/2) with the operating points.
% op: struct with rows J1, J2, phi1, phi2, alpha, xi_dw, alpha_nac1,
%     alpha_nac2 (results.dbg of a trajectory segment, or the results of
%     cruise_postp_2Wings_dbg for a steady point)
% models: models_config output; name: output PDF in Images/
% t: time of each point (optional), colours the points; omitted for a
%    single steady point

if nargin < 4 || isempty(t), t = []; end

J1 = op.J1;  J2 = op.J2;
phi1 = rad2deg(op.phi1);  phi2 = rad2deg(op.phi2);
alpha      = rad2deg(op.alpha);
alpha_tail = rad2deg(op.alpha - op.xi_dw);
alpha_nac  = rad2deg([op.alpha_nac1; op.alpha_nac2]);

if ~isempty(t) && t(end) > 120
    t = t/60; t_label = '$t$ [min]';
else
    t_label = '$t$ [s]';
end

fig = figure('Name', ['Model maps: ' name], 'Position', [100 100 1500 800]);

%% Propeller maps

J_grid   = linspace(0, 1.2, 121);
phi_grid = linspace(-90, 90, 181);
[JJ, PP] = meshgrid(J_grid, phi_grid);
maps   = {models.CT_lookup, models.CP_lookup};
titles = {'$C_T(J, \phi)$', '$C_P(J, \phi)$'};
labels = {'$C_T$', '$C_P$'};
for m = 1:2
    subplot(2,4,(m-1)*2 + (1:2)); hold on; box on;
    Z = eval_map(maps{m}, JJ, deg2rad(PP));
    levels = nice_levels(Z, 12);
    % min(Z) as first level so the region below the lowest round level is filled
    [C, hc] = contourf(JJ, PP, Z, [min(Z(:)), levels], 'LineColor', [0.55 0.55 0.55], 'HandleVisibility', 'off');
    clabel(C, hc, levels, 'FontSize', 9, 'Color', [0.3 0.3 0.3], 'LabelSpacing', 500);
    contour(JJ, PP, Z, [0 0], 'k', 'LineWidth', 2, 'DisplayName', [labels{m} ' = 0']);
    colormap(gca, 0.5 + 0.5*gray(256));            % Light grey: the coloured markers are the trajectory
    colorbar;
    plot(J1, phi1, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(J2, phi2, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, 'HandleVisibility', 'off');
    plot_points(J1, phi1, 'o', t, 'rotor 1');
    plot_points(J2, phi2, 's', t, 'rotor 2');
    if numel(J1) > 1
        text(J1(1), phi1(1), '  $t_0$', 'FontSize', 11);
        text(J1(end), phi1(end), '  $t_f$', 'FontSize', 11);
    end
    xlabel('$J$'); ylabel('$\phi$ [deg]'); title(titles{m});
    xlim([0 1.2]); ylim([-90 90]);
    legend('Location', 'northeast');
end

%% Aerodynamic polars

surfaces = {
    'Wing',      models.CL_wing_lookup,     models.CD_wing_lookup,     alpha,      '$\alpha$';
    'Tail',      models.CL_tail_lookup,     models.CD_tail_lookup,     alpha_tail, '$\alpha - \xi_{dw}$';
    'Fuselage',  models.CL_fuselage_lookup, models.CD_fuselage_lookup, alpha,      '$\alpha$';
    'Nacelles',  [],                        models.CD_nac_lookup,      alpha_nac,  '$\alpha_{nac}$'};

for s = 1:size(surfaces, 1)
    subplot(2,4,4 + s); hold on; grid on; box on;
    a_pts = surfaces{s,4};
    a_lim = [min(-30, min(a_pts(:)) - 10), max(45, max(a_pts(:)) + 10)];
    a_lim = max(min(a_lim, 180), -180);
    a_grid = linspace(a_lim(1), a_lim(2), 300);
    if ~isempty(surfaces{s,2})
        CL = eval_map(surfaces{s,2}, deg2rad(a_grid));
        plot(a_grid, CL, 'Color', [0 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', '$C_L$');
        plot_points(a_pts, eval_map(surfaces{s,2}, deg2rad(a_pts)), 'o', t, '');
    end
    CD = eval_map(surfaces{s,3}, deg2rad(a_grid));
    plot(a_grid, CD, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', '$C_D$');
    if s < 4
        plot_points(a_pts, eval_map(surfaces{s,3}, deg2rad(a_pts)), 's', t, '');
    else
        plot_points(a_pts(1,:), eval_map(surfaces{s,3}, deg2rad(a_pts(1,:))), 'o', t, 'nacelle 1');
        plot_points(a_pts(2,:), eval_map(surfaces{s,3}, deg2rad(a_pts(2,:))), 's', t, 'nacelle 2');
    end
    xlim(a_lim);
    xlabel([surfaces{s,5} ' [deg]']); title(surfaces{s,1});
    legend('Location', 'best');
end

if ~isempty(t)
    annotation('textbox', [0.3 0.005 0.4 0.03], 'String', ['Markers coloured by ' t_label ', from dark blue (start) to yellow (end)'], ...
               'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
end

picturewidth = 60;
hw_ratio = 0.5;
set(findall(fig,'-property','FontSize'),'FontSize',14)
set(findall(fig,'-property','Interpreter'),'Interpreter','latex')
set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(fig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
targetFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Images');
if ~isfolder(targetFolder), mkdir(targetFolder); end
print(fig, fullfile(targetFolder, name),'-dpdf','-vector','-fillpage')

end

function levels = nice_levels(Z, n)
% About n contour levels at round values (1, 2 or 5 x 10^k steps)
span = max(Z(:)) - min(Z(:));
raw  = span/n;
mag  = 10^floor(log10(raw));
steps = [1 2 5 10]*mag;
step  = steps(find(steps >= raw, 1));
levels = (ceil(min(Z(:))/step):floor(max(Z(:))/step))*step;
end

function Z = eval_map(f, varargin)
% Evaluate a scalar CasADi model element-wise on equally sized arrays
sz   = size(varargin{1});
args = cellfun(@(a) reshape(a, 1, []), varargin, 'UniformOutput', false);
fm   = f.map(prod(sz));
Z    = reshape(full(fm(args{:})), sz);
end

function plot_points(x, y, marker, t, label)
% Operating points, coloured by time (parula, dark blue = start, yellow =
% end, as RGB so the axes colormap of the maps is not used) or plain black
% for a single point / no time
vis = 'on';
if isempty(label), vis = 'off'; end
if isempty(t) || numel(x) == 1
    plot(x, y, marker, 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', ...
         'LineStyle', 'none', 'DisplayName', label, 'HandleVisibility', vis);
else
    cmap = parula(256);
    k    = 1 + round(255*(t - t(1))/max(t(end) - t(1), eps));
    scatter(x, y, 55, cmap(k,:), marker, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
            'DisplayName', label, 'HandleVisibility', vis);
end
end
