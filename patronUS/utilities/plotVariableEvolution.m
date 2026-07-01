function plotVariableEvolution(X_hist, bounds)
    iterations = 0:(size(X_hist,2)-1);
    idx = struct('V',1,'gamma',2,'alpha',3,'epsilon',4,'n',5);
    labels = struct('V','V [m/s]','gamma','\gamma [deg]','alpha','\alpha [deg]', ...
                     'epsilon','\epsilon [deg]','n','n [RPS]');
    is_angle = struct('V',false,'gamma',true,'alpha',true,'epsilon',true,'n',false);
    fields = fieldnames(idx);

    figure('Name','Variable Evolution: Initial to Optimal','Position',[100 100 1100 950]);
    for k = 1:numel(fields)
        f = fields{k};
        i = idx.(f);
        vals = X_hist(i,:);
        if is_angle.(f)
            vals = rad2deg(vals);
            lb = rad2deg(bounds.(f).min);
            ub = rad2deg(bounds.(f).max);
        else
            lb = bounds.(f).min;
            ub = bounds.(f).max;
        end

        subplot(5,1,k); hold on; grid on;
        plot(iterations, vals, 'b-', 'LineWidth', 1.8);
        plot(iterations(1), vals(1), 'ko', 'MarkerFaceColor','y', 'MarkerSize',8);   % initial
        plot(iterations(end), vals(end), 'ko', 'MarkerFaceColor','g', 'MarkerSize',8); % optimal
        yline(lb, 'r--', 'LineWidth', 1.2);
        yline(ub, 'r--', 'LineWidth', 1.2);
        title([f ' evolution'], 'Interpreter','none');
        ylabel(labels.(f));
        if k == numel(fields), xlabel('IPOPT Iteration'); end
    end
    legend('trajectory','initial guess','optimal','bounds','Location','bestoutside');
end