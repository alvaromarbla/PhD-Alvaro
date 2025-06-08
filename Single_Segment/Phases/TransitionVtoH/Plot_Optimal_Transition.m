function Plot_Optimal_Transition(opt_slope, base_params)
    % Run simulation with optimal parameters
    optimized_params = base_params;
    optimized_params.epsilon_slope = opt_slope;
    [results, ~] = Transition_Solver(optimized_params);
    
    % Extract derivative results 
    results.L = CL_Model(base_params.alpha_opt_eff) * 0.5*base_params.rho*results.V.^2*base_params.S_ref;
    results.D = CD_Model(base_params.alpha_opt_eff) * 0.5*base_params.rho*results.V.^2*base_params.S_ref;

    [results.T, ~] = CT_Model(results.V, results.n, base_params.alpha_opt_eff, results.epsilon_vec, base_params);
    results.P = CP_Model(results.V, results.n, base_params.alpha_opt_eff, results.epsilon_vec, base_params)*base_params.rho*base_params.prop.num_engines.*results.n.^3*base_params.prop.diameter^5;

    results.alpha_vec = ones(length(results.t))*base_params.alpha_opt_eff;
   
    % Create figure with subplots
    figure('Name', 'Optimal Transition Profile', 'Position', [100 100 1200 800]);
    hold on 
    
    % Plot 1: Tilt Angle and Velocity
    subplot(3,2,1);
    yyaxis left
    plot(results.t, rad2deg(results.epsilon_vec), 'b', 'LineWidth', 1.5);
    ylabel('Engine Tilt (deg)');
    yyaxis right
    plot(results.t, results.V, 'r', 'LineWidth', 1.5);
    hold on
    yline(base_params.V_min_ope, '--k')
    hold off
    ylabel('Velocity (m/s)');
    title(sprintf('Optimal Tilt Rate: %.1f deg/s', rad2deg(opt_slope)));
    grid on
    legend('Tilt Angle', 'Velocity', 'Location', 'northwest');
    
    % Plot 2: Energy Consumption
    subplot(3,2,2);
    plot(results.t, results.E/base_params.E_batt*100, 'LineWidth', 1.5);
    title('Energy Consumption Profile');
    xlabel('Time (s)');
    ylabel('Percentage of energy stored (%)');
    grid on
    
    % Plot 3: Engine RPM
    subplot(3,2,3);
    plot(results.t, results.n, 'LineWidth', 1.5);
    title('Engine Speed Profile');
    xlabel('Time (s)');
    ylabel('RPS');
    yline(base_params.prop.rps_max, '--r', 'Max RPM');
    grid on
    
    % Plot 4: Aerodynamic Forces
    subplot(3,2,4);
    plot(results.t, results.L, 'b', 'LineWidth', 1.5);
    hold on
    plot(results.t, results.D, 'r', 'LineWidth', 1.5);
    title('Aerodynamic Forces');
    xlabel('Time (s)');
    ylabel('Force (N)');
    legend('Lift', 'Drag', 'Location', 'northwest');
    grid on
    
    % Plot 5: Angle of Attack
    subplot(3,2,5);
    plot(results.t, rad2deg(results.alpha_vec), 'LineWidth', 1.5);
    title('Angle of Attack Profile');
    xlabel('Time (s)');
    ylabel('AoA (deg)');
    grid on
    
    % Plot 6: Thrust and Power
    subplot(3,2,6);
    yyaxis left
    plot(results.t, results.T, 'b', 'LineWidth', 1.5);
    ylabel('Thrust (N)');
    yyaxis right
    plot(results.t, results.P/1000, 'r', 'LineWidth', 1.5);
    title('Thrust and Power');
    xlabel('Time (s)');
    ylabel('Power (kW)');
    legend('Thrust', 'Power', 'Location', 'northwest');
    grid on
end