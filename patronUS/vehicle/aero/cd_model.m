function CD_lookup = cd_model(CD_version)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% DRAG COEFFICIENT MODEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version: XXX.X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import casadi.*

grid_alpha = linspace(deg2rad(-35), deg2rad(50), 100);

if CD_version == "bspline"
    % 9th-order polynomial for C_D(a) from wind tunnel data
    p_CD_ac = [0.08034; 0.0179749; 3.715949; -0.852857; -2.376796; ...
        1.049514; 0.53177; -0.485949; -0.0586846; 0.0778734];

    % Manual polynomial evaluation preserving original order
    CD_data = p_CD_ac(10)*grid_alpha.^9 + p_CD_ac(9)*grid_alpha.^8 + p_CD_ac(8)*grid_alpha.^7 + ...
        p_CD_ac(7)*grid_alpha.^6 + p_CD_ac(6)*grid_alpha.^5 + p_CD_ac(5)*grid_alpha.^4 + ...
        p_CD_ac(4)*grid_alpha.^3 + p_CD_ac(3)*grid_alpha.^2 + p_CD_ac(2)*grid_alpha + p_CD_ac(1);

    CD_lookup = interpolant('CD_lookup','bspline', {grid_alpha}, CD_data);

elseif CD_version == "analytical"

    alpha = SX.sym('alpha');
    p_CD_ac = [0.08034; 0.0179749; 3.715949; -0.852857; -2.376796; ...
        1.049514; 0.53177; -0.485949; -0.0586846; 0.0778734];

    CD_expr = p_CD_ac(10)*alpha^9 + p_CD_ac(9)*alpha^8 + p_CD_ac(8)*alpha^7 + ...
        p_CD_ac(7)*alpha^6 + p_CD_ac(6)*alpha^5 + p_CD_ac(5)*alpha^4 + ...
        p_CD_ac(4)*alpha^3 + p_CD_ac(3)*alpha^2 + p_CD_ac(2)*alpha + p_CD_ac(1);

    % Compile into a clean CasADi Function object
    CD_lookup = Function('CD_lookup', {alpha}, {CD_expr});

else
    warning('CD version unknown')

end