function CM_fuselage_lookup = cm_fuselage(CM_version)
import casadi.*
if CM_version == "analytical"
    alpha_mesh = deg2rad([-180, -170, -160, -150, -140, -130, -120, -110, -100, -90, ...
        -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, ...
        20, 30, 40, 50, 60, 70, 80, 90, 100, 110, ...
        120, 130, 140, 150, 160, 170, 180]);

    CM_fuselage_2wing_mesh = [-0.0233, 0.0256, 0.105, 0.228, 0.356, 0.474, 0.560, 0.545, 0.404, 0.383, ...
    0.371, 0.358, 0.410, 0.382, 0.327, 0.248, 0.151, 0.0680, 0.0161, 0.0147, ...
    -0.00704, -0.0971, -0.214, -0.355, -0.491, -0.448, -0.428, -0.456, -0.483, -0.506, ...
    -0.587, -0.696, -0.607, -0.465, -0.286, -0.0941, -0.0233];

    CM_fuselage_interp = interpolant('CM_fuselage_interp', 'bspline', {alpha_mesh}, CM_fuselage_2wing_mesh);

    CM_fuselage_lookup = CM_fuselage_interp;

elseif CM_version == "dummy"

    alpha = SX.sym('alpha');
    CM0 = 0.0; CM_alpha = 0.10;
    CM_expr = CM0 + CM_alpha*alpha;
    CM_fuselage_lookup = Function('CM_fuselage_lookup', {alpha}, {CM_expr});
else
    warning('CL version unknown')
end
end