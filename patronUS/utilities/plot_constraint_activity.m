function plot_constraint_activity(activity, max_rows)
% Constraint activity per node (activity from report_constraint_activity).
% One panel per constraint with at least one active node, sorted by
% sum(lambda): left axis, the constrained quantity and its limit (bounds in
% physical units, path constraints as the normalised c <= 0); right axis,
% the multiplier at each node. A node with lambda ~ 0 has slack; where
% lambda > 0 the constraint binds, and its size says how much relaxing the
% limit there would lower the objective.

if nargin < 2
    max_rows = 8;
end

sel = activity([activity.n_active] > 0);
sel = sel(1:min(end, max_rows));
if isempty(sel)
    fprintf('plot_constraint_activity: no active constraints.\n');
    return
end
n = numel(sel);

consActFig = figure('Name', 'Constraint activity per node', 'Position', [100 100 1100 230*n]);

for k = 1:n
    r = sel(k);
    subplot(n, 1, k); hold on; grid on;

    yyaxis left
    plot(r.t, r.value, '-o', 'LineWidth', 1.5, 'MarkerSize', 3);
    yline(r.limit, '--', 'LineWidth', 1.2);
    if strcmp(r.kind, 'bound')
        ylabel(to_latex(strtok(r.name)));
    else
        ylabel('$c$');
    end

    yyaxis right
    stem(r.t, r.lam, 'filled', 'MarkerSize', 3);
    ylabel('$\lambda$');

    title(sprintf('%s: $\\Sigma\\lambda$ = %.3g, active at %d/%d nodes', ...
                  to_latex(r.name), r.lam_sum, r.n_active, r.n_points));
end
xlabel('$t$ [s]');

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
print(consActFig, fullfile(targetFolder, 'constraint_activity'),'-dpdf','-vector','-fillpage')

end

function s = to_latex(s)
% Constraint names (e.g. 'T1 <= T_max') as LaTeX text
s = strrep(s, '_', '\_');
s = strrep(s, '<=', '$\leq$');
s = strrep(s, '>=', '$\geq$');
end
