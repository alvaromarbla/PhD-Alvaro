function PlotForcesAngles(sol, params)
    % Extract parameters
    V = sol(1);
    gamma = sol(2);  % Flight path angle (airflow direction)
    alpha = sol(3);  % Angle of attack (vehicle vs airflow)
    epsilon = sol(4); % Thrust vector angle (engine vs vehicle)
    n = params.n_ope;
    [T, phi] = CT_Model(V, n, alpha, epsilon, params);
    L = CL_Model(alpha) * 0.5 * params.rho * V^2 * params.wing_area;
    D = CD_Model(alpha) * 0.5 * params.rho * V^2 * params.wing_area;
    W = params.mass * params.g;
    
    % Create figure
    figure('Name', 'Force/Angle Diagram', 'Position', [100 100 800 800]);
    axis equal;
    hold on;
    
    % Reference point (vehicle location)
    origin = [0 0];
    
    %% 1. Draw coordinate system
    plot([-1.5 1.5], [0 0], 'k--', 'LineWidth', 0.5); % Horizontal
    plot([0 0], [-1.5 1.5], 'k--', 'LineWidth', 0.5); % Vertical
    text(1.4, 0.1, 'Horiz', 'HorizontalAlignment', 'right');
    text(0.1, 1.4, 'Vert', 'VerticalAlignment', 'top');
    
    %% 2. Draw forces (scaled for visualization)
    max_force = max([L, D, W, T]);
    scale = 1/max_force;
    
    % Weight (always downward)
    quiver(0, 0, 0, -W*scale, 'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', 'k');
    
    % Airflow direction (gamma)
    airflow_dir = [cos(gamma), sin(gamma)];  % Modified for descent
    plot([0 airflow_dir(1)], [0 airflow_dir(2)], 'b--', 'LineWidth', 1.5);
    
    % Lift (perpendicular to airflow)
    lift_dir = [-airflow_dir(2), airflow_dir(1)]; % Rotate 90° CCW
    quiver(0, 0, lift_dir(1)*L*scale, lift_dir(2)*L*scale, ...
          'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', 'g');
    
    % Drag (opposite airflow)
    drag_dir = -airflow_dir;
    quiver(0, 0, drag_dir(1)*D*scale, drag_dir(2)*D*scale, ...
          'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', 'r');
    
    % Thrust direction (vehicle axis + epsilon)
    vehicle_axis = airflow_dir + [lift_dir(1)*alpha, lift_dir(2)*alpha]; % Approximate
    thrust_dir = [cos(phi), sin(phi)];  % phi = alpha + epsilon
    quiver(0, 0, thrust_dir(1)*T*scale, thrust_dir(2)*T*scale, ...
          'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', 'm');
    
    %% 3. Draw angles
    angle_scale = 0.5;
    
    % Gamma (airflow to horizontal)
    DrawAngleArc(0, gamma, angle_scale, 'b', '\gamma');
    
    % Alpha (vehicle to airflow)
    DrawAngleArc(gamma, alpha, angle_scale*0.8, 'r', '\alpha');
    
    % Epsilon (thrust to vehicle)
    DrawAngleArc(gamma + alpha, epsilon, angle_scale*0.6, 'm', '\epsilon');
    
    % Phi (total thrust angle)
    text(thrust_dir(1)*T*scale*0.5, thrust_dir(2)*T*scale*0.5, ...
        sprintf('\\phi = %.1f°', rad2deg(phi)), 'Color', 'm');
    
    %% 4. Labels and legend
    text(-0.1, -W*scale-0.1, 'Weight', 'Color', 'k', 'HorizontalAlignment', 'center');
    text(lift_dir(1)*L*scale, lift_dir(2)*L*scale, ' Lift', 'Color', 'g');
    text(drag_dir(1)*D*scale, drag_dir(2)*D*scale, ' Drag', 'Color', 'r');
    text(thrust_dir(1)*T*scale, thrust_dir(2)*T*scale, ' Thrust', 'Color', 'm');
    title(sprintf('Force Diagram @ V=%.1fm/s\n', V));
    axis off;
    
    %% Nested function for angle drawing
    function DrawAngleArc(start_ang, delta_ang, radius, color, label)
        theta = linspace(start_ang, start_ang + delta_ang, 50);
        [x, y] = pol2cart(theta, radius);
        plot(x, y, 'Color', color, 'LineWidth', 1.5);
        [tx, ty] = pol2cart(start_ang + delta_ang/2, radius*1.1);
        text(tx, ty, label, 'Color', color, 'FontSize', 12);
    end
end