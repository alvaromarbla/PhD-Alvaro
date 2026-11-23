function plot_constraint_activity(activity, max_rows, name)
% Constraint activity per node (activity from report_constraint_activity).
% One panel per constraint with at least one active node, sorted by
% sum(lambda): left axis, the constrained quantity and its limit (bounds and
% actuator rates in physical units, path constraints as the normalised
% c <= 0); right axis,
% the multiplier at each node. A node with lambda ~ 0 has slack; where
% lambda > 0 the constraint binds, and its size says how much relaxing the
% limit there would lower the objective.
% With no active constraint, the rows closest to their limits are drawn.
% max_rows: most panels (default 8); name: output PDF in Images/
% (default 'constraint_activity'). Time in minutes for segments > 2 min.

if nargin < 2 || isempty(max_rows), max_rows = 8; end
if nargin < 3, name = 'constraint_activity'; end

sel = activity([activity.n_active] > 0);
none_active = isempty(sel);
if none_active
    % Interior optimum (e.g. a cruise set by the trim equalities): show the
    % path limits and actuator rates closest to their limits instead, to
    % see the margins (state bounds next to a boundary condition, such as
    % x >= 0 one node after x(0) = 0, would crowd them out)
    fprintf('plot_constraint_activity: no active constraints; plotting the path limits closest to their limits.\n');
    sel = activity(ismember({activity.kind}, {'path', 'rate'}));
    [~, order] = sort([sel.slack_min]);
    sel = sel(order);
end
sel = sel(1:min(end, max_rows));
if isempty(sel), return; end
n = numel(sel);

t_end = max(arrayfun(@(r) max(r.t), sel));
t_div = 1; t_label = '$t$ [s]';
if t_end > 120, t_div = 60; t_label = '$t$ [min]'; end

consActFig = figure('Name', ['Constraint activity: ' name], 'Position', [100 100 1100 230*n]);

for k = 1:n
    r = sel(k);
    subplot(n, 1, k); hold on; grid on;

    yyaxis left
    plot(r.t/t_div, r.value, '-o', 'LineWidth', 1.5, 'MarkerSize', 3);
    yline(r.limit, '--', 'LineWidth', 1.2);
    if any(strcmp(r.kind, {'bound', 'rate'}))
        ylabel(to_latex(strtok(r.name)));
    else
        ylabel('$c$');
    end

    yyaxis right
    stem(r.t/t_div, r.lam, 'filled', 'MarkerSize', 3);
    ylabel('$\lambda$');

    title(sprintf('%s: $\\Sigma\\lambda$ = %.3g, active at %d/%d nodes, min slack %.3g', ...
                  to_latex(r.name), r.lam_sum, r.n_active, r.n_points, r.slack_min));
end
if none_active
    sgtitle('No active inequality: constraints closest to their limits', 'Interpreter', 'latex');
end
xlabel(t_label);

picturewidth = 50;
hw_ratio = 0.2*n;
set(findall(consActFig,'-property','FontSize'),'FontSize',16)
set(findall(consActFig,'-property','Interpreter'),'Interpreter','latex')
set(findall(consActFig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(consActFig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(consActFig,'Position');
set(consActFig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
targetFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Images');
if ~isfolder(targetFolder), mkdir(targetFolder); end
print(consActFig, fullfile(targetFolder, name),'-dpdf','-vector','-fillpage')

end

function s = to_latex(s)
% Constraint names (e.g. 'T1 <= T_max') as LaTeX text
s = strrep(s, '_', '\_');
s = strrep(s, '<=', '$\leq$');
s = strrep(s, '>=', '$\geq$');
end
