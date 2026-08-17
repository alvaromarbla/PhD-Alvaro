function plotforceangles(ax, results)
    % Extract parameters
    V = results.opt.V;
    gamma = results.opt.gamma;
    alpha = results.opt.alpha;
    epsilon = results.opt.epsilon;
    alpha_nac = results.alpha_nac;

    L = results.L_opt;
    D = results.D_opt;
    T = results.T_opt;
    W = results.W_opt;
    D_nac = results.D_nac_opt;
    phi = alpha + epsilon;

    axes(ax); % Focus on this subplot axis for standard plotting behaviors
    axis(ax, 'equal');
    hold(ax, 'on');

    % Reference point
    plot(ax, [-1.5 1.5], [0 0], 'k--', 'LineWidth', 0.5);
    plot(ax, [0 0], [-1.5 1.5], 'k--', 'LineWidth', 0.5);
    text(ax, 1.4, 0.1, 'Horiz');
    text(ax, 0.1, 1.4, 'Vert');

    % Force Scaling
    max_force = max([L, D, W, T]);
    scale = 1 / max_force;

    % Draw Quivers and Lines
    quiver(ax, 0, 0, 0, -W*scale, 'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', 'k');

    airflow_dir = [cos(gamma), sin(gamma)];
    plot(ax, [0 airflow_dir(1)], [0 airflow_dir(2)], 'b--', 'LineWidth', 1.5);

    lift_dir = [-airflow_dir(2), airflow_dir(1)];
    quiver(ax, 0, 0, lift_dir(1)*L*scale, lift_dir(2)*L*scale, 'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', 'g');

    drag_dir = -airflow_dir;
    quiver(ax, 0, 0, drag_dir(1)*D*scale, drag_dir(2)*D*scale, 'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', 'r');

    drag_nac_dir = [cos(phi+pi-alpha_nac), sin(phi+pi-alpha_nac)];
    quiver(0, 0, drag_nac_dir(1)*D_nac*scale, drag_nac_dir(2)*D_nac*scale, 'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', 'r');

    thrust_dir = [cos(phi), sin(phi)];
    quiver(ax, 0, 0, thrust_dir(1)*T*scale, thrust_dir(2)*T*scale, 'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', 'm');

    % Local Arc Drawing
    angle_scale = 0.5;
    DrawLocalArc(ax, 0, gamma, angle_scale, 'b', '\gamma');
    DrawLocalArc(ax, gamma, alpha, angle_scale*0.8, 'r', '\alpha');
    DrawLocalArc(ax, gamma + alpha, epsilon, angle_scale*0.6, 'm', '\epsilon');

    text(ax, thrust_dir(1)*T*scale*0.5, thrust_dir(2)*T*scale*0.5, sprintf('\\phi = %.1f°', rad2deg(phi)), 'Color', 'm');
    text(ax, -0.1, -W*scale-0.1, 'Weight', 'Color', 'k', 'HorizontalAlignment', 'center');
    text(ax, lift_dir(1)*L*scale, lift_dir(2)*L*scale, ' Lift', 'Color', 'g');
    text(ax, drag_dir(1)*D*scale, drag_dir(2)*D*scale, ' Drag', 'Color', 'r');
    text(ax, drag_nac_dir(1)*D_nac*scale, drag_nac_dir(2)*D_nac*scale, ' Drag_{nac}', 'Color', 'r');
    text(ax, thrust_dir(1)*T*scale, thrust_dir(2)*T*scale, ' Thrust', 'Color', 'm');

    title(ax, sprintf('Force Diagram @ V=%.1fm/s', V));
    axis(ax, 'off');
end

function DrawLocalArc(ax, start_ang, delta_ang, radius, color, label)
    theta = linspace(start_ang, start_ang + delta_ang, 50);
    [x, y] = pol2cart(theta, radius);
    plot(ax, x, y, 'Color', color, 'LineWidth', 1.5);
    [tx, ty] = pol2cart(start_ang + delta_ang/2, radius*1.1);
    text(ax, tx, ty, label, 'Color', color, 'FontSize', 12);
end