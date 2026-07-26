function plotIpoptconvergence(history_matrix, bounds)
    if isempty(history_matrix)
        error('No history matrix data found. Ensure the solver executed successfully.');
    end

    iterations = 0:(size(history_matrix, 2) - 1);
    fields = fieldnames(bounds.vars);        % row k <-> fields{k}, same as X build
    n = numel(fields);

    % Which variables to show in degrees, and axis labels
    deg_vars = {'gamma', 'alpha', 'epsilon'};

    figure('Name', 'IPOPT Boundary Proximity Trace', 'Position', [100 100 1100 200*n]);

    for k = 1:n
        f  = fields{k};
        lb = bounds.vars.(f).min;
        ub = bounds.vars.(f).max;
        y  = history_matrix(k, :);

        use_deg = ismember(f, deg_vars);
        if use_deg
            y  = rad2deg(y);
            lb = rad2deg(lb);
            ub = rad2deg(ub);
        end

        subplot(n, 1, k); hold on; grid on;
        plot(iterations, y, 'b-', 'LineWidth', 2);
        plot(iterations, ones(size(iterations)) * lb, 'r--', 'LineWidth', 1.5);
        plot(iterations, ones(size(iterations)) * ub, 'r--', 'LineWidth', 1.5);
        title(sprintf('%s Evolution vs. Bounds', f), 'Interpreter', 'none');
        ylabel(f, 'Interpreter', 'none');
    end

    xlabel('IPOPT Iteration Number');
end