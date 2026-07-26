function plotforcesangles_imag(ax, results)

    % Extract parameters from your optimization pipeline
    V = results.V_opt;
    gamma = results.gamma_opt;
    alpha = results.alpha_opt;
    epsilon = results.epsilon_opt;
    alpha_nac = results.alpha_nac;

    L = results.L_opt;
    D = results.D_opt;
    T = results.T_opt;
    W = results.W_opt;
    D_nac = results.D_nac_opt;
    phi = alpha + epsilon;

    % Core Axis Setup
    axes(ax);
    cla(ax); % Clear previous iteration data for interactive updates
    axis(ax, 'equal');
    hold(ax, 'on');

    %% 1. IMAGE IMPORT & ROTATION (HGTRANSFORMS)
    % Replace these with your actual transparent CAD PNG filenames
    % Images must be cropped tightly to their boundaries
    try
        [img_fuse, ~, alpha_fuse] = imread('fuselage_sideview.png');
        [img_nac, ~, alpha_nac_img] = imread('nacelle_sideview.png');
    catch
        % Fallback placeholders for development testing
        img_fuse = uint8(200 * ones(100, 200, 3)); alpha_fuse = ones(100, 200);
        img_nac = uint8(150 * ones(40, 60, 3)); alpha_nac_img = ones(40, 60);
    end

    % Physical scaling boundaries for your CAD images
    fuse_width = 2.0; fuse_height = 1.0;
    nac_width = 0.6;  nac_height = 0.3;

    % Define wing attachments relative to fuselage CAD center (0,0)
    % Example configuration for a front and rear tandem wing layout
    wing_offsets = [
        0.4,  0.1;  % Wing 1 (Front Pair) [X, Y]
        -0.5, 0.15  % Wing 2 (Rear Pair) [X, Y]
    ];

    % Parent Fuselage Transform (Rotates with aircraft Alpha)
    t_fuse = hgtransform('Parent', ax);
    img_h1 = imagesc(ax, [-fuse_width/2, fuse_width/2], [fuse_height/2, -fuse_height/2], img_fuse, 'Parent', t_fuse);
    set(img_h1, 'AlphaData', alpha_fuse);

    % Apply overall pitch angle rotation (alpha) to fuselage
    set(t_fuse, 'Matrix', makehgtform('zrotate', alpha));

    % Child Nacelle Transforms (Attach to Fuselage, rotate by alpha_nac)
    for i = 1:2
        % Create the transform handle directly as a child of t_fuse
        t_nac(i) = hgtransform('Parent', t_fuse);

        img_n = imagesc(ax, [-nac_width/2, nac_width/2], [nac_height/2, -nac_height/2], img_nac, 'Parent', t_nac(i));
        set(img_n, 'AlphaData', alpha_nac_img);

        % Translate to wing mount, then apply local thrust tilt angle
        set(t_nac(i), 'Matrix', makehgtform('translate', [wing_offsets(i,:), 0]) * makehgtform('zrotate', alpha_nac));
    end

    %% 2. FORCE VECTOR PLOTTING
    % Reference Grid Lines
    plot(ax, [-2 2], [0 0], 'k--', 'LineWidth', 0.5);
    plot(ax, [0 0], [-2 2], 'k--', 'LineWidth', 0.5);
    text(ax, 1.8, 0.1, 'Horiz', 'FontSize', 9);

    % Scaling pipeline
    max_force = max([L, D, W, T]);
    scale = 0.8 / max_force;

    % Airframe Global forces (Origin 0,0)
    quiver(ax, 0, 0, 0, -W*scale, 'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', 'k');

    airflow_dir = [cos(gamma), sin(gamma)];
    plot(ax, [0 airflow_dir(1)], [0 airflow_dir(2)], 'b--', 'LineWidth', 1.2);

    lift_dir = [-airflow_dir(2), airflow_dir(1)];
    quiver(ax, 0, 0, lift_dir(1)*L*scale, lift_dir(2)*L*scale, 'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', [0 0.6 0]);

    drag_dir = -airflow_dir;
    quiver(ax, 0, 0, drag_dir(1)*D*scale, drag_dir(2)*D*scale, 'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', 'r');

    % Loop over the 2 wing stations to plot thrust and local nacelle drag
    for i = 1:2
        % Calculate global coordinate position of the rotated wing tip
        % This ensures arrows track perfectly with the image rotation
        R_matrix = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
        global_wing_pos = R_matrix * wing_offsets(i,:)';
        wx = global_wing_pos(1);
        wy = global_wing_pos(2);

        % Local Thrust Vector (Acting at wing mount location)
        thrust_dir = [cos(phi), sin(phi)];
        quiver(ax, wx, wy, thrust_dir(1)*(T/2)*scale, thrust_dir(2)*(T/2)*scale, 'LineWidth', 2, 'MaxHeadSize', 0.3, 'Color', 'm');

        % Nacelle Drag Vector
        drag_nac_dir = [cos(phi+pi-alpha_nac), sin(phi+pi-alpha_nac)];
        quiver(ax, wx, wy, drag_nac_dir(1)*(D_nac/2)*scale, drag_nac_dir(2)*(D_nac/2)*scale, ...
            'LineWidth', 1.5, 'MaxHeadSize', 0.3, 'Color', [0.8 0.3 0]);
    end

    %% 3. ARCS & LABELS
    angle_scale = 0.4;
    DrawLocalArc(ax, 0, gamma, angle_scale, 'b', '\gamma');
    DrawLocalArc(ax, gamma, alpha, angle_scale*0.8, 'r', '\alpha');
    DrawLocalArc(ax, gamma + alpha, epsilon, angle_scale*0.6, 'm', '\epsilon');

    % Text Callouts
    text(ax, 0, -W*scale-0.15, 'Weight', 'Color', 'k', 'HorizontalAlignment', 'center');
    text(ax, lift_dir(1)*L*scale, lift_dir(2)*L*scale, ' Total Lift', 'Color', [0 0.6 0]);

    title(ax, sprintf('Interactive Aerodynamic State: V = %.1f m/s', V));
    xlim(ax, [-2.2, 2.2]);
    ylim(ax, [-2.2, 2.2]);
    axis(ax, 'on');
    grid(ax, 'on');
end

function DrawLocalArc(ax, start_ang, delta_ang, radius, color, label)
    theta = linspace(start_ang, start_ang + delta_ang, 50);
    [x, y] = pol2cart(theta, radius);
    plot(ax, x, y, 'Color', color, 'LineWidth', 1.5);
    [tx, ty] = pol2cart(start_ang + delta_ang/2, radius*1.2);
    text(ax, tx, ty, label, 'Color', color, 'FontSize', 11, 'FontWeight', 'bold');
end