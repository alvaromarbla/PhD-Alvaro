function generate_results_page_2Wings(results, models,params)
    % Create a series of master figures
    figure('Name', 'Patronus Performance Propulsion Results Page', 'Position', [100, 100, 1200, 900]);

    J_opt1 = results.J1;
    phi_opt1 =  results.phi1;
    J_opt2 = results.J2;
    phi_opt2 =  results.phi2;

    alpha_opt = results.opt.alpha;
    xi_dw     = results.xi_dw;

    CL_wing_opt = results.CL_wing;
    CD_wing_opt = results.CD_wing;
    CM_wing_opt = results.CM_wing;

    CL_fuselage_opt = results.CL_fuselage;
    CD_fuselage_opt = results.CD_fuselage;
    CM_fuselage_opt = results.CM_fuselage;

    CL_tail_opt = results.CL_tail;
    CD_tail_opt = results.CD_tail;
    CM_tail_opt = results.CM_tail;

    alpha_nac1 = results.alpha_nac1;
    CD_nac1_opt = results.CD_nac1;
    alpha_nac2 = results.alpha_nac2;
    CD_nac2_opt = results.CD_nac2;

    V = results.opt.V;
    vi1 = results.vi1;
    vi2 = results.vi2;
    V_nac1 = results.V_nac1;
    V_nac2 = results.V_nac2;
    phi1 = results.phi1;
    phi2 = results.phi2;
    %% ====================================================
    %% SUBPLOT 1: CP1 Contour Map
    %% ====================================================
    ax1 = subplot(3, 2, 1);
    plot_cp_contour(ax1, J_opt1,phi_opt1, models);

    %% ====================================================
    %% SUBPLOT 2: CP2 Contour Map
    %% ====================================================
    ax2 = subplot(3, 2, 2);
    plot_cp_contour(ax2, J_opt2,phi_opt2, models);

    %% ====================================================
    %% SUBPLOT 3: CT1 Contour Map
    %% ====================================================
    ax3 = subplot(3, 2, 3);
    plot_ct_contour(ax3, J_opt1,phi_opt1, models);

    %% ====================================================
    %% SUBPLOT 4: CT2 Contour Map
    %% ====================================================
    ax4 = subplot(3, 2, 4);
    plot_ct_contour(ax4, J_opt2,phi_opt2, models);

    %% ====================================================
    %% SUBPLOT 5: CH1 Contour Map
    %% ====================================================
    ax5 = subplot(3, 2, 5);
    plot_ch_contour(ax5, J_opt1,phi_opt1, models);

    %% ====================================================
    %% SUBPLOT 6: CH2 Contour Map
    %% ====================================================
    ax6 = subplot(3, 2, 6);
    plot_ch_contour(ax6, J_opt2,phi_opt2, models);

    figure('Name', 'Patronus Performance Aero Wing Results Page', 'Position', [100, 100, 1200, 900]);

    %% ====================================================
    %% SUBPLOT 1: Aerodynamic Curves (CL vs Alpha)
    %% ====================================================
    ax1 = subplot(2, 2, 1);
    plot_CL_wing(ax1, alpha_opt, CL_wing_opt, models);

    %% ====================================================
    %% SUBPLOT 2: Aerodynamic Curves (CD vs Alpha)
    %% ====================================================
    ax2 = subplot(2, 2, 2);
    plot_CD_wing(ax2, alpha_opt, CD_wing_opt, models);

    %% ====================================================
    %% SUBPLOT 3: Aerodynamic Curves (CD vs Alpha)
    %% ====================================================
    ax3 = subplot(2, 2, 3);
    plot_CM_wing(ax3, alpha_opt, CM_wing_opt, models);

    %% ====================================================
    %% SUBPLOT 4: Aerodynamic Curves (CL vs Alpha)
    %% ====================================================
    ax4 = subplot(2, 2, 4);
    plot_AeroEff_wing(ax4, alpha_opt, CL_wing_opt, CD_wing_opt, models);

    figure('Name', 'Patronus Performance Aero Fuselage Results Page', 'Position', [100, 100, 1200, 900]);

    %% ====================================================
    %% SUBPLOT 1: Aerodynamic Curves (CL vs Alpha)
    %% ====================================================
    ax1 = subplot(2, 2, 1);
    plot_CL_fuselage(ax1, alpha_opt, CL_fuselage_opt, models);

    %% ====================================================
    %% SUBPLOT 2: Aerodynamic Curves (CD vs Alpha)
    %% ====================================================
    ax2 = subplot(2, 2, 2);
    plot_CD_fuselage(ax2, alpha_opt, CD_fuselage_opt, models);

     %% ====================================================
    %% SUBPLOT 3: Aerodynamic Curves (CD vs Alpha)
    %% ====================================================
    ax3 = subplot(2, 2, 3);
    plot_CM_fuselage(ax3, alpha_opt, CM_fuselage_opt, models);

    %% ====================================================
    %% SUBPLOT 4: Aerodynamic Curves (Aero Eff)
    %% ====================================================
    ax4 = subplot(2, 2, 4);
    plot_AeroEff_fuselage(ax4, alpha_opt, CL_fuselage_opt, CD_fuselage_opt, models);

    figure('Name', 'Patronus Performance Aero Tail Results Page', 'Position', [100, 100, 1200, 900]);

    %% ====================================================
    %% SUBPLOT 1: Aerodynamic Curves (CL vs Alpha)
    %% ====================================================
    ax1 = subplot(2, 2, 1);
    plot_CL_tail(ax1, alpha_opt-xi_dw, CL_tail_opt, models);

    %% ====================================================
    %% SUBPLOT 2: Aerodynamic Curves (CD vs Alpha)
    %% ====================================================
    ax2 = subplot(2, 2, 2);
    plot_CD_tail(ax2,  alpha_opt-xi_dw, CD_tail_opt, models);

     %% ====================================================
    %% SUBPLOT 3: Aerodynamic Curves (CD vs Alpha)
    %% ====================================================
    ax3 = subplot(2, 2, 3);
    plot_CM_tail(ax3,  alpha_opt-xi_dw, CM_tail_opt, models);

    %% ====================================================
    %% SUBPLOT 4: Aerodynamic Curves (Aero Eff)
    %% ====================================================
    ax4 = subplot(2, 2, 4);
    plot_AeroEff_tail(ax4,  alpha_opt-xi_dw, CL_tail_opt, CD_tail_opt, models);


    %% ====================================================
    %% PLOT 5: Nacelle speed triangles and drag
    %% ====================================================
    figure('Name', 'Patronus Nacelles Results Page', 'Position', [100, 100, 1200, 900]);

    %% SUBPLOT 1 and 2: Speed Triangle
    axTriangle1 = subplot(2, 2, 1);
    plotspeedtriangle_nac(axTriangle1, V , vi1, V_nac1, phi1);

    axTriangle2 = subplot(2, 2, 2);
    plotspeedtriangle_nac(axTriangle2, V , vi2, V_nac2, phi2);

    %% SUBPLOT 3 and 4: CD for Nacelle

    ax3 = subplot(2, 2, 3);
    plot_CD_nac(ax3, alpha_nac1, CD_nac1_opt, models);

    ax4 = subplot(2, 2, 4);
    plot_CD_nac(ax4, alpha_nac2, CD_nac2_opt, models);

    %% ====================================================
    %% PLOT 6: Forces and Angles Diagram
    %% ====================================================

    figure('Name', 'Patronus Forces Results Page', 'Position', [100, 100, 1200, 900]);

    axForce = axes();
    %plotforceangles(axForce, results); % This works for 1Wing because I dont have the picture yet. To be improved.
    plotforcesangles_imag(axForce, results,params);
end

%% ================================================================
%% HELPER FUNCTIONS (Appended at the bottom of the same file)
%% ================================================================

function plot_cp_contour(ax, J_opt,phi_opt, models)
    % Extract parameters
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

function plot_ct_contour(ax, J_opt,phi_opt , models)
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
    hold(ax, "on");

    % Highlight 0-line
    contour(ax, J_vec, phi_vec*180/pi, CT_mat', 'LineWidth', 2, 'LevelList', 0, 'LineColor', "r");

    % Highlight Operating Point
    plot(ax, J_opt, phi_opt*180/pi, 'ro', 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'LineWidth', 2);

    grid(ax, 'on');
    ylabel(ax, 'Engine Tilt Angle $\phi$ [deg]', 'Interpreter', 'latex');
    xlabel(ax, 'Advance Ratio J [-]');
    title(ax, '$C_T$ Model Map \& Optimal Point', 'Interpreter', 'latex');
end

function plot_ch_contour(ax, J_opt,phi_opt, models)
    % Extract parameters
    CH = models.CH_lookup;

    J_vec = linspace(0, 1.2, 50);
    phi_vec = linspace(0, pi/2, 90);
    CH_mat = zeros(length(J_vec), length(phi_vec));

    for ii = 1:length(J_vec)
        for jj = 1:length(phi_vec)
            CH_mat(ii,jj) = full(CH(J_vec(ii), phi_vec(jj)));
        end
    end

    % Plot contour onto specified axis
    [c, h] = contourf(ax, J_vec, phi_vec*180/pi, CH_mat', ...
        'LevelList', [-0.05 -0.04 -0.03 -0.02 -0.01 0 0.01 0.02 0.026 0.028 0.03 0.035 0.04]);
    clabel(c, h, "Interpreter", "latex");
    hold(ax, 'on');

    % Highlight 0-line and Operating Point
    contour(ax, J_vec, phi_vec*180/pi, CH_mat', 'LineWidth', 2, 'LevelList', 0, 'LineColor', "r");
    plot(ax, J_opt, phi_opt*180/pi, 'ro', 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'LineWidth', 2);

    grid(ax, 'on');
    ylabel(ax, 'Engine Tilt Angle $\phi$ [deg]', 'Interpreter', 'latex');
    xlabel(ax, 'Advance Ratio J [-]');
    title(ax, '$C_H$ Model Map \& Optimal Point', 'Interpreter', 'latex');
end

function plot_CL_wing(ax, alpha_opt, CL_opt, models)

    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)
    % Example linear lift curve: CL = CL0 + CLa * alpha
    CL_vec = full(models.CL_wing_lookup(alpha_vec));
    % Plot CL
    plot(ax, rad2deg(alpha_vec), CL_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CL_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$C_L$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep for Wing: $C_L$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_CL_fuselage(ax, alpha_opt, CL_opt, models)

    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)
    % Example linear lift curve: CL = CL0 + CLa * alpha
    CL_vec = full(models.CL_fuselage_lookup(alpha_vec));
    % Plot CL
    plot(ax, rad2deg(alpha_vec), CL_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CL_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$C_L$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep for Fuselage: $C_L$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_CL_tail(ax, alpha_opt, CL_opt, models)

    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)
    % Example linear lift curve: CL = CL0 + CLa * alpha
    CL_vec = full(models.CL_tail_lookup(alpha_vec));
    % Plot CL
    plot(ax, rad2deg(alpha_vec), CL_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CL_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$C_L$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep for Tail: $C_L$ vs. $\alpha - \xi_{dw}$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_CD_wing(ax, alpha_opt, CD_opt, models)
    % Extract current optimal points
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)
    % Example linear lift curve: CL = CL0 + CLa * alpha
    CD_vec = full(models.CD_wing_lookup(alpha_vec));

    % Plot CD
    plot(ax, rad2deg(alpha_vec), CD_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CD_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$C_D$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep for Wing: $C_D$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_CD_fuselage(ax, alpha_opt, CD_opt, models)
    % Extract current optimal points
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)
    % Example linear lift curve: CL = CL0 + CLa * alpha
    CD_vec = full(models.CD_fuselage_lookup(alpha_vec));

    % Plot CD
    plot(ax, rad2deg(alpha_vec), CD_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CD_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$C_D$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep for Fuselage: $C_D$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');
end

function plot_CD_tail(ax, alpha_opt, CD_opt, models)
    % Extract current optimal points
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)
    % Example linear lift curve: CL = CL0 + CLa * alpha
    CD_vec = full(models.CD_tail_lookup(alpha_vec));

    % Plot CD
    plot(ax, rad2deg(alpha_vec), CD_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CD_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$C_D$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep for Tail: $C_D$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_CM_wing(ax, alpha_opt, CM_opt, models)
    % Extract current optimal points
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)
    % Example linear lift curve: CL = CL0 + CLa * alpha
    CM_vec = full(models.CM_wing_lookup(alpha_vec));

    % Plot CD
    plot(ax, rad2deg(alpha_vec), CM_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CM_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$C_M$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep for Wing: $C_M$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_CM_fuselage(ax, alpha_opt, CM_opt, models)
    % Extract current optimal points
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)
    % Example linear lift curve: CL = CL0 + CLa * alpha
    CM_vec = full(models.CM_fuselage_lookup(alpha_vec));

    % Plot CD
    plot(ax, rad2deg(alpha_vec), CM_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CM_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$C_M$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep for Fuselage: $C_M$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_CM_tail(ax, alpha_opt, CM_opt, models)
    % Extract current optimal points
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)
    % Example linear lift curve: CL = CL0 + CLa * alpha
    CM_vec = full(models.CM_tail_lookup(alpha_vec));

    % Plot CD
    plot(ax, rad2deg(alpha_vec), CM_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CM_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$C_M$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep for Tail: $C_M$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_AeroEff_wing(ax, alpha_opt, CL_opt,CD_opt , models)
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)

    CL_vec = full(models.CL_wing_lookup(alpha_vec));
    CD_vec = full(models.CD_wing_lookup(alpha_vec));

    % Plot Aero Eff
    plot(ax, rad2deg(alpha_vec), CL_vec./CD_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CL_opt/CD_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$Aero_{eff}$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep for Wing: $Aero_{eff}$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_AeroEff_fuselage(ax, alpha_opt, CL_opt,CD_opt , models)
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)

    CL_vec = full(models.CL_fuselage_lookup(alpha_vec));
    CD_vec = full(models.CD_fuselage_lookup(alpha_vec));

    % Plot Aero Eff
    plot(ax, rad2deg(alpha_vec), CL_vec./CD_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CL_opt/CD_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$Aero_{eff}$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep for Fuselage: $Aero_{eff}$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_AeroEff_tail(ax, alpha_opt, CL_opt,CD_opt , models)
    % Vector for full curve generation
    alpha_vec = linspace(-pi/2, pi/2, 50);

    % Dynamic CL curve calculation (Replace with your model's formula/lookup)

    CL_vec = full(models.CL_tail_lookup(alpha_vec));
    CD_vec = full(models.CD_tail_lookup(alpha_vec));

    % Plot Aero Eff
    plot(ax, rad2deg(alpha_vec), CL_vec./CD_vec, 'b-', 'LineWidth', 2);
    hold(ax, 'on');

    % Plot Optimal Operating Point (Convert alpha to degrees if necessary)
    plot(ax, rad2deg(alpha_opt), CL_opt/CD_opt, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    grid(ax, 'on');
    xlabel(ax, 'Angle of Attack $\alpha$ [deg]', 'Interpreter', 'latex');
    ylabel(ax, '$Aero_{eff}$ [-]', 'Interpreter', 'latex');
    title(ax, 'Aerodynamic Sweep for Tail: $Aero_{eff}$ vs. $\alpha$', 'Interpreter', 'latex');
    legend(ax, 'Model Curve', 'Optimal Point', 'Location', 'best');

end

function plot_CD_nac(ax, alpha_nac, CD_nac_opt, models)
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

function plotspeedtriangle_nac(ax, V_inf , v_i, V_nac, phi)

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