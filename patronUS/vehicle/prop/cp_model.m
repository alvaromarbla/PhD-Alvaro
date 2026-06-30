function CP_lookup = cp_model(CP_version)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% POWER COEFFICIENT MODEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version: XXX.X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import casadi.*

if CP_version == "bspline"

    % 1. Define regular grids for your independent variables
    % Choose limits that safely cover your bounds.V and bounds.alpha/epsilon
    grid_J   = linspace(0, 1.2, 50);          % 60 points for Advance Ratio (J)
    grid_phi = linspace(-pi/2, pi/2, 100);    % 60 points for Thrust Angle (phi in radians)

    % 2. Evaluate your original CP polynomial on this 2D grid to create the carpet data
    [Mesh_J, Mesh_Phi] = meshgrid(grid_J, grid_phi);
    CP_data = zeros(size(Mesh_J));

    % Re-use your original polynomial logic numerically (pure double math here)
    a1_P = 0.0261518307541734;   b1_P = 0.0473735972985378;   c1_P = -0.16267474946046;
    d1_P = 0.0247028469343899;   e1_P = 0.0306053713439883;   b2_P = -0.0762350484603968;
    c2_P = 0.148580471912353;    d2_P = -0.0726017200715775;  b3_P = 0.0897273366920878;
    c3_P = 0.0122602815262456;   b4_P = -0.0486029866039398;

    for i = 1:numel(Mesh_J)
        j_val = Mesh_J(i);
        p_val = Mesh_Phi(i);

        CP_data(i) = a1_P + b1_P*j_val + c1_P*j_val^2 + d1_P*j_val^3 + e1_P*j_val^4 + ...
            abs(p_val)*(b2_P*j_val + c2_P*j_val^2 + d2_P*j_val^3) + ...
            abs(p_val)^2*(b3_P*j_val + c3_P*j_val^2) + ...
            abs(p_val)^3*b4_P*j_val;
    end

    % 3. Flatten the data column-wise to match CasADi's storage requirement
    % CasADi expects a 1D vector of outputs ordered by flattening the grid
    CP_flat = CP_data(:);

    % 4. Instantiate the CasADi B-spline Function object
    % Knots are automatically generated linearly by passing degrees [3, 3] (cubic spline)
    % Name of function, grid inputs, flattened output data, [degree_dim1, degree_dim2]
    CP_lookup = interpolant('CP_lookup','bspline', {grid_J, grid_phi}, CP_flat);

elseif CP_version == "analytical"
    J = SX.sym('J');
    phi = SX.sym('phi');

    % Re-use your original polynomial logic numerically (pure double math here)
    a1_P = 0.0261518307541734;   b1_P = 0.0473735972985378;   c1_P = -0.16267474946046;
    d1_P = 0.0247028469343899;   e1_P = 0.0306053713439883;   b2_P = -0.0762350484603968;
    c2_P = 0.148580471912353;    d2_P = -0.0726017200715775;  b3_P = 0.0897273366920878;
    c3_P = 0.0122602815262456;   b4_P = -0.0486029866039398;

    abs_phi_smooth = sqrt(phi^2 + 1e-6);

    CP_expr = a1_P + b1_P*J + c1_P*J^2 + d1_P*J^3 + e1_P*J^4 + ...
        abs_phi_smooth*(b2_P*J + c2_P*J^2 + d2_P*J^3) + ...
        (abs_phi_smooth^2)*(b3_P*J + c3_P*J^2) + ...
        (abs_phi_smooth^3)*b4_P*J;

    CP_lookup = Function('CP_lookup', {J,phi}, {CP_expr});

else
    warning ('CP version unknown')

end

end