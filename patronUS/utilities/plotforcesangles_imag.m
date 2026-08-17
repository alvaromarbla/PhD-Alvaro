function plotforcesangles_imag(ax, results,params)

    % Extract parameters from your optimization pipeline
    V = results.opt.V;
    gamma = results.opt.gamma;
    alpha = results.opt.alpha;
    epsilon1 = results.opt.epsilon1;
    epsilon2 = results.opt.epsilon2;
    xi_dw  = results.xi_dw;

    L_wing = results.L_wing;
    D_wing = results.D_wing;
    L_fuselage = results.L_fuselage;
    D_fuselage = results.D_fuselage;
    L_tail = results.L_tail;
    D_tail = results.D_tail;
    alpha_nac1 = results.alpha_nac1;
    alpha_nac2 = results.alpha_nac2;

    T1 = results.T1;
    T2 = results.T2;
    W = results.W;
    D_nac1 = results.D_nac1;
    D_nac2 = results.D_nac2;
    phi1 = alpha + epsilon1;
    phi2 = alpha + epsilon2;

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
    fuse_width = 2.2; fuse_height = 0.8;
    nac_width = 0.5;  nac_height = 0.2;

    % Parent Fuselage Transformation Object (Rotates dynamically with pitch alpha)
    t_fuse = hgtransform('Parent', ax);
    img_h1 = imagesc(ax, [-fuse_width/2, fuse_width/2], [fuse_height/2, -fuse_height/2], img_fuse, 'Parent', t_fuse);
    set(img_h1, 'AlphaData', alpha_fuse);
    set(t_fuse, 'Matrix', makehgtform('zrotate', alpha));

    % Child Nacelle Transformations (Attached to physical structure positions)
    % Nacelle 1 (Wing location mount)
    t_nac1 = hgtransform('Parent', t_fuse);
    img_n1 = imagesc(ax, [-nac_width/2, nac_width/2], [nac_height/2, -nac_height/2], img_nac, 'Parent', t_nac1);
    set(img_n1, 'AlphaData', alpha_nac_img);
    set(t_nac1, 'Matrix', makehgtform('translate', [params.geo.xw, params.geo.zw, 0]) * makehgtform('zrotate', epsilon1));

    % Nacelle 2 (Tail location mount)
    t_nac2 = hgtransform('Parent', t_fuse);
    img_n2 = imagesc(ax, [-nac_width/2, nac_width/2], [nac_height/2, -nac_height/2], img_nac, 'Parent', t_nac2);
    set(img_n2, 'AlphaData', alpha_nac_img);
    set(t_nac2, 'Matrix', makehgtform('translate', [params.geo.xtw, params.geo.ztw, 0]) * makehgtform('zrotate', epsilon2));

    %% 3. Station Origin Tracking Matrices (Global Coordinates)
    % Rotational matrix tracking the body coordinate frame transformation
    R_body = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];

    % Transform CAD offset boundaries into current global frame positions
    pos_CG   = [0; 0]; % Center of mass chosen as reference origin (0,0)
    pos_wing = R_body * [params.geo.xw;  params.geo.zw];
    pos_fus  = R_body * [params.geo.xfus; params.geo.zfus];
    pos_tail = R_body * [params.geo.xtw;  params.geo.ztw];

    %% 4. Vector Force Computations & Plot Scaling Pipeline
    % Determine consistent visual layout vector scaling limits
    max_force = max([L_wing, D_wing, W, T1, T2, L_tail]);
    scale = 0.6 / max_force;

    % Baseline Freestream Vectors (Freestream airflow acts relative to gamma angle)
    dir_freestream = [cos(gamma), sin(gamma)];
    dir_free_lift  = [-dir_freestream(2), dir_freestream(1)];
    dir_free_drag  = -dir_freestream;

    % --- STATION 1: CENTER OF MASS ---
    quiver(ax, pos_CG(1), pos_CG(2), 0, -W*scale, 'LineWidth', 2.5, 'MaxHeadSize', 0.4, 'Color', 'k');
    text(ax, pos_CG(1), pos_CG(2) - W*scale - 0.1, 'Weight', 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', "center");

    % --- STATION 2: MAIN WING ---
    % Thrust 1 Vector (oriented at phi1 relative to horizontal)
    dir_T1 = [cos(phi1), sin(phi1)];
    quiver(ax, pos_wing(1), pos_wing(2), dir_T1(1)*T1*scale, dir_T1(2)*T1*scale, 'LineWidth', 2, 'Color', 'm');

    % Wing Aerodynamic Forces (align with freestream airflow)
    quiver(ax, pos_wing(1), pos_wing(2), dir_free_lift(1)*L_wing*scale, dir_free_lift(2)*L_wing*scale, 'LineWidth', 2, 'Color', [0 0.5 0]);
    quiver(ax, pos_wing(1), pos_wing(2), dir_free_drag(1)*D_wing*scale, dir_free_drag(2)*D_wing*scale, 'LineWidth', 2, 'Color', 'r');

    % Nacelle 1 Structural Drag (acts at alpha_nac1 relative to freestream local vector)
    dir_D_nac1 = -[cos(alpha_nac1), sin(alpha_nac1)];
    quiver(ax, pos_wing(1), pos_wing(2), dir_D_nac1(1)*D_nac1*scale, dir_D_nac1(2)*D_nac1*scale, 'LineWidth', 1.5, 'Color', [0.8 0.4 0]);

    % --- STATION 3: FUSELAGE AERO CENTER ---
    % Fuselage Aerodynamic Forces (align with freestream airflow)
    quiver(ax, pos_fus(1), pos_fus(2), dir_free_lift(1)*L_fuselage*scale, dir_free_lift(2)*L_fuselage*scale, 'LineWidth', 2, 'Color', [0 0.7 0]);
    quiver(ax, pos_fus(1), pos_fus(2), dir_free_drag(1)*D_fuselage*scale, dir_free_drag(2)*D_fuselage*scale, 'LineWidth', 2, 'Color', [0.9 0.1 0]);

    % --- STATION 4: TAIL STATION ---
    % Tail Local Airflow Angle (Freestream modified directly by downwash angle xi_dw)
    gamma_tail = gamma - xi_dw;
    dir_tail_lift = [-sin(gamma_tail), cos(gamma_tail)];
    dir_tail_drag = -[cos(gamma_tail), sin(gamma_tail)];

    % Tail Aerodynamic Forces (Aligned with the downwash modified airflow)
    quiver(ax, pos_tail(1), pos_tail(2), dir_tail_lift(1)*L_tail*scale, dir_tail_lift(2)*L_tail*scale, 'LineWidth', 2, 'Color', [0.1 0.6 0.3]);
    quiver(ax, pos_tail(1), pos_tail(2), dir_tail_drag(1)*D_tail*scale, dir_tail_drag(2)*D_tail*scale, 'LineWidth', 2, 'Color', [0.7 0.2 0.2]);

    % Thrust 2 Vector (oriented at phi2 relative to horizontal)
    dir_T2 = [cos(phi2), sin(phi2)];
    quiver(ax, pos_tail(1), pos_tail(2), dir_T2(1)*T2*scale, dir_T2(2)*T2*scale, 'LineWidth', 2, 'Color', [0.6 0 0.6]);

    % Nacelle 2 Structural Drag (acts at alpha_nac2 relative to local tail vector)
    dir_D_nac2 = -[cos(alpha_nac2), sin(alpha_nac2)];
    quiver(ax, pos_tail(1), pos_tail(2), dir_D_nac2(1)*D_nac2*scale, dir_D_nac2(2)*D_nac2*scale, 'LineWidth', 1.5, 'Color', [0.8 0.5 0.1]);

    %% 5. DIAGNOSTIC INTERACTION ANGLE ARCS & LABELS
    angle_scale = 0.5;
    plot(ax, [pos_CG(1), pos_CG(1)+2.0], [pos_CG(2), pos_CG(2)], 'k:', 'LineWidth', 0.8); % Horizon indicator line

    DrawLocalArc(ax, pos_CG(1), pos_CG(2), 0, gamma, angle_scale, 'b', '\gamma');
    DrawLocalArc(ax, pos_CG(1), pos_CG(2), gamma, alpha, angle_scale*0.8, 'r', '\alpha');
    DrawLocalArc(ax, pos_wing(1), pos_wing(2), alpha, epsilon1, angle_scale*0.6, 'm', '\epsilon_1');
    DrawLocalArc(ax, pos_tail(1), pos_tail(2), alpha, epsilon2, angle_scale*0.6, [0.6 0 0.6], '\epsilon_2');

    % Text callouts for identification
    text(ax, pos_wing(1), pos_wing(2)+0.15, 'Wing', 'FontSize', 8, 'FontWeight', 'bold');
    text(ax, pos_tail(1), pos_tail(2)+0.15, 'Tail', 'FontSize', 8, 'FontWeight', 'bold');

    title(ax, sprintf('Generalized Aircraft Force Distribution Profile | V = %.1f m/s', V));
    xlim(ax, [-2.5, 2.5]);
    ylim(ax, [-2.5, 2.5]);
    grid(ax, 'on');
    box(ax, 'on');
end

function DrawLocalArc(ax, ox, oy, start_ang, delta_ang, radius, color, label)
    theta = linspace(start_ang, start_ang + delta_ang, 40);
    x = ox + radius * cos(theta);
    y = oy + radius * sin(theta);
    plot(ax, x, y, 'Color', color, 'LineWidth', 1.2);

    tx = ox + radius * 1.25 * cos(start_ang + delta_ang/2);
    ty = oy + radius * 1.25 * sin(start_ang + delta_ang/2);
    text(ax, tx, ty, label, 'Color', color, 'FontSize', 10, 'FontWeight', 'bold');
end