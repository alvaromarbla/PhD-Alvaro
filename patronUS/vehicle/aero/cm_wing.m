function CM_wing_lookup = cm_wing(CM_version)
import casadi.*
if CM_version == "analytical"
    alpha_mesh = deg2rad([-180, -170, -160, -150, -140, -130, -120, -110, -100, -90, ...
        -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, ...
        20, 30, 40, 50, 60, 70, 80, 90, 100, 110, ...
        120, 130, 140, 150, 160, 170, 180]);

    CM_wing_mesh = [0.0683, 0.318, 0.325, 0.391, 0.442, 0.480, 0.508, 0.525, 0.527, 0.516, ...
    0.489, 0.451, 0.399, 0.330, 0.245, 0.141, 0.0462, -0.0606, 0.0954, -0.0541, ...
    -0.178, -0.270, -0.353, -0.434, -0.480, -0.532, -0.571, -0.599, -0.612, -0.613, ...
    -0.602, -0.581, -0.537, -0.458, -0.398, -0.364, 0.0683];

    CM_wing_interp = interpolant('CM_wing_interp', 'bspline', {alpha_mesh}, CM_wing_mesh);

    CM_wing_lookup = CM_wing_interp;

elseif CM_version == "dummy"
    alpha = SX.sym('alpha');
    CM0 = -0.05; CM_alpha = -0.25;
    CM_expr = CM0 + CM_alpha*alpha;
    CM_wing_lookup = Function('CM_wing_lookup', {alpha}, {CM_expr});

else
    warning('CL version unknown')
end
end