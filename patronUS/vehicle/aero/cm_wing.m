function CM_wing_lookup = cm_wing(CM_version)
import casadi.*
if CM_version == "analytical"
    alpha_mesh = deg2rad([-180, -170, -160, -150, -140, -130, -120, -110, -100, -90, ...
        -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, ...
        20, 30, 40, 50, 60, 70, 80, 90, 100, 110, ...
        120, 130, 140, 150, 160, 170, 180]);

    CM_wing_mesh = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ...
        ];

    CM_wing_interp = interpolant('CM_wing_interp', 'bspline', {alpha_mesh}, CM_wing_mesh);

    CM_wing_lookup = CM_wing_interp;

else
    warning('CL version unknown')
end
end