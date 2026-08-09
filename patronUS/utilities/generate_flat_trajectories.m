function [splines] = generate_flat_trajectories(V_cruise, tf_guess)
    % V_cruise: Target forward speed at end of transition (e.g., 35 m/s)
    % tf_guess: Guess for total transition duration (e.g., 20 seconds)

    %% 1. Non-dimensional Time Grid
    num_points = 100;
    tau = linspace(0, 1, num_points);

    % Define the number of nodes (control intervals) to display on the spatial plot
    num_nodes = 20;
    tau_nodes = linspace(0, 1, num_nodes);

    %% 2. Construct x(tau) Trajectory (Minimum 3rd order polynomial)
    A_x = [3, 2; 6, 2];
    B_x = [V_cruise * tf_guess; 0];
    coeffs_x = A_x \ B_x;
    poly_x = [coeffs_x(1), coeffs_x(2), 0, 0];

    %% 3. Construct z(tau) Trajectory (Minimum 5th order polynomial)
    A_z = [1,  1,  1;
           5,  4,  3;
           20, 12, 6];
    B_z = [10; 0; 0];
    coeffs_z = A_z \ B_z;
    poly_z = [coeffs_z(1), coeffs_z(2), coeffs_z(3), 0, 0, 0];

    %% 4. Construct theta(tau) Trajectory (Minimum 3rd order polynomial)
    theta_f = 3 * pi / 180;
    A_th = [1, 1; 3, 2];
    B_th = [theta_f; 0];
    coeffs_th = A_th \ B_th;
    poly_theta = [coeffs_th(1), coeffs_th(2), 0, 0];

    %% 5. Evaluate Trajectories and Derivatives
    splines.tau = tau;
    splines.t   = tau * tf_guess;

    % Positions
    splines.x     = polyval(poly_x, tau);
    splines.z     = polyval(poly_z, tau);
    splines.theta = polyval(poly_theta, tau);

    % Node positions (for visualization anchors)
    splines.x_nodes = polyval(poly_x, tau_nodes);
    splines.z_nodes = polyval(poly_z, tau_nodes);

    % First Derivatives (dx/dtau)
    poly_dx = polyder(poly_x);
    poly_dz = polyder(poly_z);
    poly_dth = polyder(poly_theta);

    % Convert to physical velocities (dx/dt = (1/tf)*dx/dtau)
    splines.vx = (1/tf_guess) * polyval(poly_dx, tau);
    splines.vz = (1/tf_guess) * polyval(poly_dz, tau);
    splines.q  = (1/tf_guess) * polyval(poly_dth, tau);

    % Second Derivatives (d2x/dtau2)
    poly_ddx = polyder(poly_dx);
    poly_ddz = polyder(poly_dz);
    poly_ddth = polyder(poly_dth);

    % Convert to physical accelerations (d2x/dt2 = (1/tf^2)*d2x/dtau2)
    splines.ax    = (1/tf_guess^2) * polyval(poly_ddx, tau);
    splines.az    = (1/tf_guess^2) * polyval(poly_ddz, tau);
    splines.q_dot = (1/tf_guess^2) * polyval(poly_ddth, tau);

    %% 6. Derive Flatness Quantities
    splines.V     = sqrt(splines.vx.^2 + splines.vz.^2);
    splines.gamma = atan2(splines.vz, splines.vx); % Flightpath angle in radians
    splines.alpha = splines.theta - splines.gamma;

    %% 7. Expanded Visualization Layout

    % --- Figure 1: Spatial & Flight Path States ---
    figure('Name', 'Spatial Trajectory and Flight Profile', 'NumberTitle', 'off');

    % Subplot 1: Downrange Distance x vs Time t
    subplot(2,2,1);
    plot(splines.t, splines.x, 'b-', 'LineWidth', 1.5);
    ylabel('Downrange Distance X (m)');
    xlabel('Time t (s)');
    grid on;
    title('Forward Position Evolution');

    % Subplot 2: Altitude Z vs Downrange X with Spline Nodes
    subplot(2,2,2);
    plot(splines.x, splines.z, 'r-', 'LineWidth', 1.5); hold on;
    plot(splines.x_nodes, splines.z_nodes, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
    ylabel('Altitude Z (m)');
    xlabel('Downrange Distance X (m)');
    legend('Continuous Path', 'Spline Nodes/Knots', 'Location', 'best');
    grid on;
    title('Spatial Flight Path Profile');

    % Subplot 3: Velocity Components (vx, vz) vs Time t
    subplot(2,2,3);
    plot(splines.t, splines.vx, 'g-', 'LineWidth', 1.5); hold on;
    plot(splines.t, splines.vz, 'm-', 'LineWidth', 1.5);
    ylabel('Velocity (m/s)');
    xlabel('Time t (s)');
    legend('V_x (Forward)', 'V_z (Vertical)', 'Location', 'best');
    grid on;
    title('Velocity Components Profile');

    % Subplot 4: Flight Path Angle (gamma) vs Time t
    subplot(2,2,4);
    plot(splines.t, splines.gamma * 180 / pi, 'k-', 'LineWidth', 1.5);
    ylabel('\gamma (degrees)');
    xlabel('Time t (s)');
    grid on;
    title('Flight Path Angle Evolution');

    % --- Figure 2: Attitude & Velocity Tracking (Supplementary) ---
    figure('Name', 'Attitude and Total Velocity Profile', 'NumberTitle', 'off');
    subplot(2,1,1);
    plot(splines.t, splines.V, 'LineWidth', 1.5);
    ylabel('Total Velocity V (m/s)');
    grid on;
    title('Total Kinematic Speed');

    subplot(2,1,2);
    plot(splines.t, splines.theta*180/pi, 'b-', 'LineWidth', 1.5); hold on;
    plot(splines.t, splines.alpha*180/pi, 'r--', 'LineWidth', 1.5);
    ylabel('Angles (deg)');
    xlabel('Time t (s)');
    legend('\theta (Pitch)', '\alpha (AoA)', 'Location', 'best');
    grid on;
    title('Attitude & Aerodynamic Angles');
end