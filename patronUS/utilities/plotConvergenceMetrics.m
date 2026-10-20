function plotConvergenceMetrics(stats, viol, slack)
    iterations = 0:(numel(stats.iterations.inf_pr)-1);

    ipoptMetFig = figure('Name','Constraint Violation & Slack','Position',[100 100 1100 800]);

    % Panel 1: IPOPT's own infeasibility metrics
    subplot(3,1,1); hold on; grid on;
    semilogy(iterations, max(stats.iterations.inf_pr, eps), 'b-', 'LineWidth', 2);
    semilogy(iterations, max(stats.iterations.inf_du, eps), 'r-', 'LineWidth', 2);
    yline(1e-9, 'k--'); % your opts.ipopt.tol
    title('IPOPT Reported Infeasibility'); ylabel('Infeasibility (log)');
    legend('inf\_pr (primal)','inf\_du (dual)','tol','Location','best');

    % Panel 2: your hand-computed constraint violation
    subplot(3,1,2); hold on; grid on;
    semilogy(iterations, max(viol.ineq_max, eps), 'm-', 'LineWidth', 1.8);
    semilogy(iterations, max(viol.eq_max, eps), 'c-', 'LineWidth', 1.8);
    title('Constraint Violation (per type)'); ylabel('Violation (log)');
    legend('max ineq violation','max eq violation','Location','best');

    % Panel 3: bound slack, normalized, per variable
    subplot(3,1,3); hold on; grid on;
    fields = fieldnames(slack);
    colors = lines(numel(fields));
    for k = 1:numel(fields)
        semilogy(iterations, max(slack.(fields{k}).normalized, eps), ...
                 'Color', colors(k,:), 'LineWidth', 1.5);
    end
    title('Normalized Bound Slack (0 = at bound)'); ylabel('Slack (log)');
    xlabel('IPOPT Iteration'); legend(fields, 'Location','best');

    picturewidth = 50;
    hw_ratio = 0.75;
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex')
    set(findall(ipoptMetFig,'-property','FontSize'),'FontSize',16)
    set(findall(ipoptMetFig,'-property','Box'),'Box','off') % optional
    set(findall(ipoptMetFig,'-property','Interpreter'),'Interpreter','latex')
    set(findall(ipoptMetFig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    set(ipoptMetFig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
    pos = get(ipoptMetFig,'Position');
    set(ipoptMetFig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    targetFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Images');
    if ~isfolder(targetFolder), mkdir(targetFolder); end
    outputpath = fullfile(targetFolder, 'ipoptMetric');
    print(ipoptMetFig, outputpath,'-dpdf','-vector','-fillpage')  
end