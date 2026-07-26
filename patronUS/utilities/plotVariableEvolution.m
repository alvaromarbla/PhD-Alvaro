function plotVariableEvolution(X_hist, bounds)
    iterations = 0:(size(X_hist,2)-1);
    fields = fieldnames(bounds.vars);       % row k <-> fields{k}
    n = numel(fields);

    % Presentation metadata, keyed by name; anything not listed prints raw
    labels   = struct('V','V [m/s]','gamma','\gamma [deg]','alpha','\alpha [deg]', ...
                      'epsilon','\epsilon [deg]','n','n [RPS]');
    is_angle = struct('gamma',true,'alpha',true,'epsilon',true);

    figure('Name','Variable Evolution: Initial to Optimal','Position',[100 100 1100 190*n]);

    for k = 1:n
        f    = fields{k};
        vals = X_hist(k,:);
        lb   = bounds.vars.(f).min;
        ub   = bounds.vars.(f).max;

        if isfield(is_angle, f) && is_angle.(f)
            vals = rad2deg(vals);
            lb   = rad2deg(lb);
            ub   = rad2deg(ub);
        end

        subplot(n,1,k); hold on; grid on;
        plot(iterations, vals, 'b-', 'LineWidth', 1.8);
        plot(iterations(1),   vals(1),   'ko', 'MarkerFaceColor','y', 'MarkerSize',8);
        plot(iterations(end), vals(end), 'ko', 'MarkerFaceColor','g', 'MarkerSize',8);
        yline(lb, 'r--', 'LineWidth', 1.2);
        yline(ub, 'r--', 'LineWidth', 1.2);
        title([f ' evolution'], 'Interpreter','none');

        if isfield(labels, f)
            ylabel(labels.(f));
        else
            ylabel(f, 'Interpreter','none');
        end

        if k == n, xlabel('IPOPT Iteration'); end
    end

    legend('trajectory','initial guess','optimal','bounds','Location','bestoutside');
end