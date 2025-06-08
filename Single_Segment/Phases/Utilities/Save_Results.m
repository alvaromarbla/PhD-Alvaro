function Save_Results(sol, params)
    % Saves optimization results in .mat and CSV formats with metadata
    % Inputs:
    %   sol    - Optimization solution vector
    %   params - UAV configuration parameters
    
    %% Create output directory if needed
    output_dir = fullfile(pwd, 'Optimization_Results');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    %% Create timestamped filename
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    base_name = sprintf('VTOL_Descent_%s', timestamp);
    mat_file = fullfile(output_dir, [base_name '.mat']);
    csv_file = fullfile(output_dir, [base_name '.csv']);
    
    %% Extract and compute key parameters
    results = struct();
    
    % Optimization variables
    results.V = sol(1);               % Airspeed [m/s]
    results.gamma = sol(2);           % Flight path angle [rad]
    results.alpha = sol(3);           % Angle of attack [rad]
    results.epsilon = sol(4);         % Thrust vector angle [rad]
    results.n = sol(5);               % Motor RPM
    
    % Derived parameters
    [results.T, results.phi] = CT_Model(results.V, results.n, ...
                                      results.alpha, results.epsilon, params);
    results.L = CL_Model(results.alpha) * 0.5 * params.rho * results.V^2 * params.wing_area;
    results.D = CD_Model(results.alpha) * 0.5 * params.rho * results.V^2 * params.wing_area;
    results.CP = CP_Model(results.V, results.n, results.alpha, results.epsilon, params);
    results.J = results.V / (results.n * params.prop.diameter);  % Advance ratio
    results.xf = -params.deltaH * cot(results.gamma)  % Horizontal distance
    
    % Force ratios
    results.L_over_W = results.L / (params.mass * params.g);
    results.T_over_W = results.T / (params.mass * params.g);
    results.L_over_D = results.L / results.D;
    
    % Metadata
    results.timestamp = timestamp;
    results.software_version = 'VTOL Optimizer 1.0';
    results.param_config = params;  % Save full parameter set
    
    %% Save in MATLAB format
    save(mat_file, '-struct', 'results');
    fprintf('Results saved to MAT file: %s\n', mat_file);
    
    %% Optional CSV Export (for external analysis)
    try
        % Create table of key parameters
        csv_data = struct2table(struct(...
            'V', results.V, ...
            'gamma_deg', rad2deg(results.gamma), ...
            'alpha_deg', rad2deg(results.alpha), ...
            'epsilon_deg', rad2deg(results.epsilon), ...
            'n_rpm', results.n, ...
            'Thrust_N', results.T, ...
            'Lift_N', results.L, ...
            'Drag_N', results.D, ...
            'CP', results.CP, ...
            'Advance_Ratio', results.J, ...
            'Horizontal_Distance_m', results.xf, ...
            'L_over_W', results.L_over_W, ...
            'T_over_W', results.T_over_W, ...
            'L_over_D', results.L_over_D ...
        ));
        
        writetable(csv_data, csv_file);
        fprintf('Results saved to CSV file: %s\n', csv_file);
    catch ME
        warning(ME.identifier, 'CSV export failed: %s', ME.message);
    end
    
    %% Create README file with format description
    readme_file = fullfile(output_dir, 'RESULTS_FORMAT.md');
    if ~exist(readme_file, 'file')
        fid = fopen(readme_file, 'w');
        fprintf(fid, '# VTOL Optimization Results Format\n\n');
        fprintf(fid, '## MAT File Contents\n');
        fprintf(fid, '- **V**: Airspeed [m/s]\n');
        fprintf(fid, '- **gamma**: Flight path angle [rad]\n');
        fprintf(fid, '- **alpha**: Angle of attack [rad]\n');
        fprintf(fid, '- **epsilon**: Thrust vector angle [rad]\n');
        fprintf(fid, '- **n**: Motor RPM\n');
        fprintf(fid, '- **T**: Total thrust [N]\n');
        fprintf(fid, '- **phi**: Effective thrust angle [rad]\n');
        fprintf(fid, '- **L**: Lift force [N]\n');
        fprintf(fid, '- **D**: Drag force [N]\n');
        fprintf(fid, '- **CP**: Power coefficient\n');
        fprintf(fid, '- **J**: Advance ratio\n');
        fprintf(fid, '- **xf**: Horizontal distance [m]\n');
        fprintf(fid, '- **L_over_W**: Lift-to-weight ratio\n');
        fprintf(fid, '- **T_over_W**: Thrust-to-weight ratio\n');
        fprintf(fid, '- **L_over_D**: Lift-to-drag ratio\n');
        fprintf(fid, '- **timestamp**: Simulation timestamp\n');
        fprintf(fid, '- **software_version**: Code version\n');
        fprintf(fid, '- **param_config**: Full parameter structure\n');
        fclose(fid);
    end
end