function CL_lookup = cl_model(CL_version)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% LIFT COEFFICIENT MODEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version: XXX.X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import casadi.*

grid_alpha = linspace(deg2rad(-35), deg2rad(50), 100);

if CL_version == "bspline"

    % 8th-order polynomial for C_L(a) from wind tunnel data
    p_CL_ac = [0.582; 3.345731; -0.642635; -2.187085; 0.713766; ...
        0.377985; -0.2946314; -0.004890; 0.043976];

    % Manual polynomial evaluation (matches original code structure)
    CL_data = p_CL_ac(9)*grid_alpha.^8 + p_CL_ac(8)*grid_alpha.^7 + p_CL_ac(7)*grid_alpha.^6 + ...
        p_CL_ac(6)*grid_alpha.^5 + p_CL_ac(5)*grid_alpha.^4 + p_CL_ac(4)*grid_alpha.^3 + ...
        p_CL_ac(3)*grid_alpha.^2 + p_CL_ac(2)*grid_alpha + p_CL_ac(1);

    CL_lookup = interpolant('CL_lookup','bspline', {grid_alpha}, CL_data);

elseif CL_version == "analytical"

    alpha = SX.sym('alpha');
    % 8th-order polynomial for C_L(a) from wind tunnel data
    p_CL_ac = [0.582; 3.345731; -0.642635; -2.187085; 0.713766; ...
        0.377985; -0.2946314; -0.004890; 0.043976];

    % Manual polynomial evaluation (matches original code structure)
    CL_expr = p_CL_ac(9)*alpha.^8 + p_CL_ac(8)*alpha.^7 + p_CL_ac(7)*alpha.^6 + ...
        p_CL_ac(6)*alpha.^5 + p_CL_ac(5)*alpha.^4 + p_CL_ac(4)*alpha.^3 + ...
        p_CL_ac(3)*alpha.^2 + p_CL_ac(2)*alpha + p_CL_ac(1);

    % Compile into a clean CasADi Function object
    CL_lookup = Function('CL_lookup', {alpha}, {CL_expr});

else
    warning ('CL version unknown')

end

end