function CT_lookup = ct_model(CT_version)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% THRUST COEFFICIENT MODEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version: XXX.X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import casadi.*

if CT_version == "bspline"

    % 1. Define regular grids for your independent variables
    % Choose limits that safely cover your bounds.V and bounds.alpha/epsilon
    grid_J   = linspace(0, 1.2, 50);          % 60 points for Advance Ratio (J)
    grid_phi = linspace(-pi/2, pi/2, 100);    % 60 points for Thrust Angle (phi in radians)

    % 2. Evaluate your original CP polynomial on this 2D grid to create the carpet data
    [Mesh_J, Mesh_Phi] = meshgrid(grid_J, grid_phi);
    CT_data = zeros(size(Mesh_J));

    % Re-use your original polynomial logic numerically (pure double math here)
    a1_T = 0.0735880531010883;   b1_T = -0.0311758018412727;   c1_T = -0.249744726429543;
    d1_T = 0.143084420694372;    e1_T = 0.0261032283758581;    b2_T = -0.0982459868751664;
    c2_T = 0.20127470719351;     d2_T = -0.173738749783189;    b3_T = 0.156239779501715;
    c3_T = 0.0368592048084175;   b4_T = -0.0478709034281346;

    for i = 1:numel(Mesh_J)
        j_val = Mesh_J(i);
        p_val = Mesh_Phi(i);

        CT_data(i) = a1_T + b1_T*j_val + c1_T*j_val^2 + d1_T*j_val^3 + e1_T*j_val^4 + ...
            abs(p_val)*(b2_T*j_val + c2_T*j_val^2 + d2_T*j_val^3) + ...
            abs(p_val)^2*(b3_T*j_val + c3_T*j_val^2) + ...
            abs(p_val)^3*b4_T*j_val;
    end

    % 3. Flatten the data column-wise to match CasADi's storage requirement
    % CasADi expects a 1D vector of outputs ordered by flattening the grid
    CT_flat = CT_data(:)';

    % 4. Instantiate the CasADi B-spline Function object
    % Knots are automatically generated linearly by passing degrees [3, 3] (cubic spline)
    % Name of function, grid inputs, flattened output data, [degree_dim1, degree_dim2]
    CT_lookup = interpolant('CT_lookup', 'bspline', {grid_J, grid_phi}, CT_flat);

elseif CT_version == "analytical"

    J = SX.sym('J');
    phi = SX.sym('phi');

    % Re-use your original polynomial logic numerically (pure double math here)
    a1_T = 0.0735880531010883;   b1_T = -0.0311758018412727;   c1_T = -0.249744726429543;
    d1_T = 0.143084420694372;    e1_T = 0.0261032283758581;    b2_T = -0.0982459868751664;
    c2_T = 0.20127470719351;     d2_T = -0.173738749783189;    b3_T = 0.156239779501715;
    c3_T = 0.0368592048084175;   b4_T = -0.0478709034281346;

    abs_phi_smooth = sqrt(phi^2 + 1e-6);

    CT_expr = a1_T + b1_T*J + c1_T*J^2 + d1_T*J^3 + e1_T*J^4 + ...
        abs_phi_smooth*(b2_T*J + c2_T*J^2 + d2_T*J^3) + ...
        (abs_phi_smooth^2)*(b3_T*J + c3_T*J^2) + ...
        (abs_phi_smooth^3)*b4_T*J;

    CT_lookup = Function('CT_lookup', {J,phi}, {CT_expr});

else
    warning ('CT version unknown')

end

end