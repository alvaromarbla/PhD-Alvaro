function generate_results_page(results, models)
    % Create a single master results figure window
    figure('Name', 'Patronus Performance Results Page', 'Position', [100, 100, 1200, 900]);

    %% ==============================================
    %% SUBPLOT 1: CP Contour Map
    %% ==============================================
    ax1 = subplot(3, 2, 1);
    plot_cp_contour(ax1, results, models);

    %% ==============================================
    %% SUBPLOT 2: CT Contour Map
    %% ==============================================
    ax2 = subplot(3, 2, 2);
    plot_ct_contour(ax2, results, models);

    %% ==============================================
    %% SUBPLOT 3: Aerodynamic Curves (CL vs Alpha)
    %% ==============================================
    ax3 = subplot(3, 2, 3);
    plot_CL(ax3, results, models);

    %% ==============================================
    %% SUBPLOT 4: Aerodynamic Curves (CL vs Alpha)
    %% ==============================================
    ax4 = subplot(3, 2, 4);
    plot_CD(ax4, results, models);

    %% ==============================================
    %% SUBPLOT 5: Aerodynamic Curves (CL vs Alpha)
    %% ==============================================
    ax5 = subplot(3, 2, 5);
    plot_AeroEff(ax5, results, models);

    %% ==============================================
    %% SUBPLOT 4: CD for Nacelle
    %% ==============================================

    ax6 = subplot(3, 2, 6);
    plot_CD_nac(ax6, results, models);

    %% ==============================================
    %% PLOT 5: Forces and Angles Diagram
    %% ==============================================

    figure('Name', 'Patronus Forces Results Page', 'Position', [100, 100, 1200, 900]);

    axForce = axes();
    plotforceangles(axForce, results);
    %plotforcesangles_imag(axForce, results);

    %% ==============================================
    %% PLOT 6: Speed triangle of nacelle
    %% ==============================================
    figure('Name', 'Patronus Speed Triangle Results Page', 'Position', [100, 100, 1200, 900]);

    axTriangle = axes();
    plotspeedtriangle_nac(axTriangle, results);
end

%% ============================================================
%% HELPER FUNCTIONS (Appended at the bottom of the same file)
%% ============================================================

function plot_cp_contour(ax, results, models)
    % Extract parameters
    J_opt = results.J_opt;
    phi_opt = results.phi_opt;
    CP = models.CP_lookup;

    J_vec = linspace(0, 1.2, 50);
    phi_vec = linspace(0, pi/2, 90);
    CP_mat = zeros(length(J_vec), length(phi_vec));

    for ii = 1:length(J_vec)
        for jj = 1:length(phi_vec)
            CP_mat(ii,jj) = full(CP(J_vec(ii), phi_vec(jj)));
        end
    end

    % Plot contour onto specified axis
    [c, h] = contourf(ax, J_vec, phi_vec*180/pi, CP_mat', ...
        'LevelList', [-0.05 -0.04 -0.03 -0.02 -0.01 0 0.01 0.02 0.026 0.028 0.03 0.035 0.04]);
    clabel(c, h, "Interpreter", "latex");
    hold(ax, 'on');

    % Highlight 0-line and Operating Point
    contour(ax, J_vec, phi_vec*180/pi, CP_mat', 'LineWidth', 2, 'LevelList', 0, 'LineColor', "r");
    plot(ax, J_opt, phi_opt*180/pi, 'ro', 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'LineWidth', 2);

    grid(ax, 'on');
    ylabel(ax, 'Engine Tilt Angle $\phi$ [deg]', 'Interpreter', 'latex');
    xlabel(ax, 'Advance Ratio J [-]');
    title(ax, '$C_P$ Model Map \& Optimal Point', 'Interpreter', 'latex');
end

function plot_ct_contour(ax, results, models)
    % Extract parameters
    J_opt = results.J_opt;
    phi_opt = results.phi_opt;
    CT = models.CT_lookup;

    J_vec = linspace(0, 1.2, 50);
    phi_vec = linspace(0, pi/2, 90);
    CT_mat = zeros(length(J_vec), length(phi_vec));

    for ii = 1:length(J_vec)
        for jj = 1:length(phi_vec)
            CT_mat(ii,jj) = full(CT(J_vec(ii), phi_vec(jj)));
        end
    end

    % Plot contour onto specified axis
    [c, h] = contourf(ax, J_vec, phi_vec*180/pi, CT_mat', 20);
    clabel(c, h, "Interpreter", "latex");
    hold(ax, 'on');

    % Highlight 0-line
    contour(ax, J_vec, phi_vec*180/pi, CT_mat', 'LineWidth', 2, 'LevelList', 0, 'LineColor', "r");

    % Highlight Operating Point
    plot(ax, J_opt, phi_opt*180/pi, 'ro', 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'LineWidth', 2);

    grid(ax, 'on');
    ylabel(ax, 'Engine Tilt Angle $\phi$ [deg]', 'Interpreter', 'latex');
    xlabel(ax, 'Advance Ratio J [-]');
    title(ax, '$C_T$ Model Map \& Optimal Point', 'Interpreter', 'latex');
end

function plot_CL(ax, results, models)
    % Extract current optimal points
    alpha_opt = results.alpha_opt; % Assuming in radians
    CL_opt = results.CL_opt;       % Optimal CL
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)
    % Example linear lift curve: CL = CL0 + CLa * alpha
    CL_vec = full(models.CL_lookup(alpha_vec));
    % Plot CL
    plot(ax, rad2deg(alpha_vec), CL_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CL_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$C_L$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep: $C_L$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_CD(ax, results, models)
    % Extract current optimal points
    alpha_opt = results.alpha_opt; % Assuming in radians
    CD_opt = results.CD_opt;
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)
    % Example linear lift curve: CL = CL0 + CLa * alpha
    CD_vec = full(models.CD_lookup(alpha_vec));

    % Plot CD
    plot(ax, rad2deg(alpha_vec), CD_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CD_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$C_D$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep: $C_D$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_AeroEff(ax, results, models)
    % Extract current optimal points
    alpha_opt = results.alpha_opt; % Assuming in radians
    CL_opt = results.CL_opt;
    CD_opt = results.CD_opt;
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)

    CL_vec = full(models.CL_lookup(alpha_vec));
    CD_vec = full(models.CD_lookup(alpha_vec));

    % Plot CD
    plot(ax, rad2deg(alpha_vec), CL_vec./CD_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CL_opt/CD_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$Aero_{eff}$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep: $Aero_{eff}$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_CD_nac(ax, results, models)
    % Extract current optimal points
    alpha_nac = results.alpha_nac; % Assuming in radians
    CD_nac_opt = results.CD_nac_opt;
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)
    % Example linear lift curve: CL = CL0 + CLa * alpha
    CD_nac_vec = full(models.CD_fuselage_lookup(alpha_vec));

    % Plot CD
    plot(ax, rad2deg(alpha_vec), CD_nac_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_nac), CD_nac_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$C_D$ Nacelle [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep: $C_D$ Nacelle vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'North');

end

function plotspeedtriangle_nac(ax, results)

    % Extract parameters
    V_inf = results.V_opt;       % Freestream velocity magnitude
    v_i   = results.vi;          % Induced velocity magnitude
    V_nac = results.V_nac;       % Total nacelle velocity magnitude
    phi   = results.phi_opt;     % Engine tilt angle (rad)
    %alpha_nac = results.alpha_nac;   % Vehicle angle of attack (rad)

    %% 1. Define Vector Components (Wind Axis System)
    % Vector A: Freestream Velocity (Horizontal)
    V_inf_vec = [-V_inf, 0];

    % Vector B: Induced Velocity (Acts down/back relative to the rotor plane)
    % Usually vi acts normal to the rotor disk. If disk is tilted by phi:
    % Normal vector pointing rearward/downward is [-sin(phi), -cos(phi)]
    % Adjust these trig components if your model defines downwash direction differently!
    vi_vec = [-v_i * sin(phi), -v_i * cos(phi)];

    % Vector C: Total Nacelle Velocity (Head-to-tail sum)
    V_nac_vec = V_inf_vec + vi_vec;

    axes(ax); % Focus on this subplot axis for standard plotting behaviors
    axis(ax, 'equal');
    hold(ax, 'on');

    %% 2. Plot the Vectors Head-to-Tail
    % Origin point for the triangle
    orig = [0,0 ] ;

    % Plot V_inf (Blue)
    quiver(ax, orig(1), orig(2), V_inf_vec(1), V_inf_vec(2), 0, ...
        'Color', [0 0.4470 0.7410], 'LineWidth', 2.5, 'MaxHeadSize', 0.2);

    % Plot v_i starting from the head of V_inf (Red)
    quiver(ax, V_inf_vec(1), V_inf_vec(2), vi_vec(1), vi_vec(2), 0, ...
        'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5, 'MaxHeadSize', 0.2);

    % Plot V_nac from origin to the final head (Green)
    quiver(ax, orig(1), orig(2), V_nac_vec(1), V_nac_vec(2), 0, ...
        'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.5, 'MaxHeadSize', 0.2);

    %% 3. Add Labels and Text Annotations
    % Calculate text placement midpoints to avoid overlap
    mid_inf = V_inf_vec / 2;
    mid_vi  = V_inf_vec + (vi_vec / 2);
    mid_nac = V_nac_vec / 2;

    text(ax, mid_inf(1), mid_inf(2) + 1.5, sprintf('V_{\\infty} = %.1f m/s', V_inf), ...
        'Color', [0 0.4470 0.7410], 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    text(ax, mid_vi(1) + 1.0, mid_vi(2), sprintf('v_i = %.1f m/s', v_i), ...
        'Color', [0.8500 0.3250 0.0980], 'FontWeight', 'bold', 'HorizontalAlignment', 'left');

    text(ax, mid_nac(1) - 1.5, mid_nac(2) - 1.0, sprintf('V_{nac} = %.1f m/s', V_nac), ...
        'Color', [0.4660 0.6740 0.1880], 'FontWeight', 'bold', 'HorizontalAlignment', 'right');

    %% 4. Clean up Title and Limits
    title(ax, 'Nacelle Velocity Inflow Triangle', 'Interpreter', 'tex');
    xlabel(ax, 'Horizontal Velocity [m/s]');
    ylabel(ax, 'Vertical Velocity [m/s]');

    % Give the plot some breathing room padding
    xlim(ax, [min([0, V_inf_vec(1), V_nac_vec(1)]) - 3, max([0, V_inf_vec(1), V_nac_vec(1)]) + 3]);
    ylim(ax, [min([0, V_inf_vec(2), vi_vec(2), V_nac_vec(2)]) - 3, max([0, V_inf_vec(2), V_nac_vec(2)]) + 3]);
end