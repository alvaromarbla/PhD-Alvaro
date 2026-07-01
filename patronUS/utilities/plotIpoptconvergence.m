function plotIpoptconvergence(history_matrix, bounds)
    % Verify the matrix actually contains iteration data
    if isempty(history_matrix)
        error('No history matrix data found. Ensure the solver executed successfully.');
    end

    % Calculate the iteration vector starting at 0
    iterations = 0:(size(history_matrix, 2) - 1);

    % Map indices based on your X state vector structure:
    % X(1)=V, X(2)=gamma, X(3)=alpha, X(4)=epsilon, X(5)=n
    idx_V   = 1;
    idx_gam = 2;
    idx_alp = 3;
    idx_eps = 4;
    idx_n   = 5;

    figure('Name', 'IPOPT Boundary & Constraint Proximity Trace', 'Position', [100 100 1100 850]);

    % --- Panel 1: Airspeed (V) ---
    subplot(4,1,1); hold on; grid on;
    V_history = history_matrix(idx_V, :);
    plot(iterations, V_history, 'b-', 'LineWidth', 2);
    % Draw bounds directly from your struct
    plot(iterations, ones(size(iterations)) * bounds.V.min, 'r--', 'LineWidth', 1.5);
    plot(iterations, ones(size(iterations)) * bounds.V.max, 'r--', 'LineWidth', 1.5);
    title('Airspeed (V) Evolution vs. Bounds');
    ylabel('V [m/s]');

    % --- Panel 2: Angle of Attack (\alpha) ---
    subplot(4,1,2); hold on; grid on;
    alpha_history = rad2deg(history_matrix(idx_alp, :));
    plot(iterations, alpha_history, 'g-', 'LineWidth', 2);
    % Draw bounds converted to degrees for clean visual interpretation
    plot(iterations, ones(size(iterations)) * rad2deg(bounds.alpha.min), 'r--', 'LineWidth', 1.5);
    plot(iterations, ones(size(iterations)) * rad2deg(bounds.alpha.max), 'r--', 'LineWidth', 1.5);
    title('Angle of Attack (\alpha) Evolution vs. Bounds');
    ylabel('\alpha [deg]');

    % --- Panel 3: Thrust Vectoring (\epsilon) ---
    subplot(4,1,3); hold on; grid on;
    eps_history = rad2deg(history_matrix(idx_eps, :));
    plot(iterations, eps_history, 'm-', 'LineWidth', 2);
    % Draw bounds converted to degrees
    plot(iterations, ones(size(iterations)) * rad2deg(bounds.epsilon.min), 'r--', 'LineWidth', 1.5);
    plot(iterations, ones(size(iterations)) * rad2deg(bounds.epsilon.max), 'r--', 'LineWidth', 1.5);
    title('Thrust Vector Angle (\epsilon) Evolution vs. Bounds');
    ylabel('\epsilon [deg]');

    % --- Panel 4: Engine Motor Speed (n) ---
    subplot(4,1,4); hold on; grid on;
    n_history = history_matrix(idx_n, :);
    plot(iterations, n_history, 'k-', 'LineWidth', 2);
    % Draw bounds directly from your struct
    plot(iterations, ones(size(iterations)) * bounds.n.min, 'r--', 'LineWidth', 1.5);
    plot(iterations, ones(size(iterations)) * bounds.n.max, 'r--', 'LineWidth', 1.5);
    title('Engine Motor Speed (n) Evolution vs. Bounds');
    xlabel('IPOPT Iteration Number');
    ylabel('n [RPS]');
end