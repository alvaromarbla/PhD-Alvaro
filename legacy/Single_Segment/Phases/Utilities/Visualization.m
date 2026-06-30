function Visualization(sol, params)
    % Extract optimized parameters
    V = sol(1);             % Airspeed [m/s]
    gamma = sol(2);         % Flight path angle [rad]
    alpha = sol(3);         % Angle of attack [rad]
    epsilon = sol(4);       % Thrust vector angle [rad]
    n = params.n_ope;             % Motor RPM
    
    % Calculate derived quantities
    [T, phi] = CT_Model(V, n, alpha, epsilon, params)
    L = CL_Model(alpha) * 0.5 * params.rho * V^2 * params.wing_area
    D = CD_Model(alpha) * 0.5 * params.rho * V^2 * params.wing_area
    CP = CP_Model(V, n, alpha, epsilon, params);
    W = params.mass * params.g;
    
    % Create figure with 2x3 subplot grid
    figure('Name', 'Optimization Results', ...
           'Position', [100 100 1200 800], ...
           'NumberTitle', 'off');
    
    %% 1. Flight Trajectory
    subplot(2,3,1);
    deltaH = params.deltaH;
    xf =  deltaH * abs(cot(gamma));
    fplot(@(h) xf - (xf/deltaH)*h, [0 deltaH], 'LineWidth', 2);
    xlabel('Altitude [m]');
    ylabel('Horizontal Distance [m]');
    title('Trajectory');
    grid on;
    axis equal tight
    
    %% 2. Force Balance Diagram
    subplot(2,3,2);
    forces = [L, T*sin(phi), D, T*cos(phi), W];
    labels = {'Lift', 'Thrust Vertical', 'Drag', 'Thrust Horizontal', 'Weight'};
    bar(forces);
    set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('Force Magnitude [N]');
    title('Force Balance');
    grid on;
    
    %% 3. Control Inputs
    subplot(2,3,3);
    controls = [...
        rad2deg(alpha), ...
        rad2deg(epsilon), ...
        rad2deg(phi), ...
        n, ...
        rad2deg(gamma)];
    
    control_labels = {...
        sprintf('AoA %.1f°', rad2deg(alpha)), ...
        sprintf('Thrust Tilt %.1f°', rad2deg(epsilon)), ...
        sprintf('Total φ  %.1f°', rad2deg(phi)), ...
        sprintf('RPS %.0f', n), ...
        sprintf('FPA %.1f°', rad2deg(gamma))};
    
    bar(controls);
    set(gca, 'XTickLabel', control_labels);
    ylabel('Control Values');
    title('Optimized Control Parameters');
    grid on;
    
    %% 4. Power Characteristics
    subplot(2,3,4);
    J = V / (n * params.prop.diameter);  % Advance ratio
    power_data = [CP, J, V/(n*D)];       % Using D as prop diameter
    bar(power_data);
    set(gca, 'XTickLabel', {'C_P', 'Advance Ratio (J)', 'V/(nD)'});
    ylabel('Dimensionless Parameters');
    title('Propulsion System Characteristics');
    grid on;
    
    %% 5. Energy Metrics
    subplot(2,3,5);
    % Energy calculation
    t_climb = params.deltaH / (V*sin(gamma))
    P = params.prop.num_engines * CP * params.rho * n^3 * params.prop.diameter^5
    Energy = params.prop.num_engines * CP * params.rho * n^3 * params.prop.diameter^5 * t_climb/params.prop.eff
    Energy_perc = Energy/ params.E_batt *100
    if params.phase == "descent"

        J_obj = (cot(gamma));
        obj_metrics = [J_obj ];     % Objective function value]
    elseif params.phase == "climb"
        x_benefit = params.deltaH * cot(gamma) * params.C_cost;
        J_obj = Energy- x_benefit;  % Minimize this
        obj_metrics = [...
        J_obj, ...      % Objective function value
        Energy, ... % Vertical force ratio
        x_benefit];    % Horizontal force ratio
    end
    
    

    
    
    bar(obj_metrics);
    set(gca, 'XTickLabel', {'J', 'Energy', 'Hor. Distance Cost'});
    ylabel('Dimensionless Ratios');
    title('Objective Function Metrics');
    grid on;
    
    %% 6. Text Summary
    subplot(2,3,6);
    axis off;
    
    summary_text = {...
        'Optimization Summary:', ...
        sprintf('Horizontal Distance: %.1f m', xf), ...
        sprintf('Flight Path Angle: %.2f°', rad2deg(gamma)), ...
        sprintf('Airspeed: %.2f m/s', V), ...
        sprintf('Thrust: %.2f N', T), ...
        sprintf('L/D Ratio: %.2f', L/D), ...
        sprintf('Power Coeff (C_P): %.2e', CP), ...
        sprintf('kJ per km @ cruise (C_{cost}): %.2e', params.C_cost), ...
        sprintf('Advance Ratio (J): %.3f', J),...
        sprintf('Engine tilt (phi): %.3f', phi*180/pi)};
    
    text(0.1, 0.8, summary_text, ...
        'VerticalAlignment', 'top', ...
        'FontSize', 10, ...
        'FontName', 'FixedWidth');

    %% 7. Forces and angles in point-mass
    PlotForcesAngles(sol, params);

end