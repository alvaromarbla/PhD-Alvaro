function CL_tail_lookup = cl_tail(CL_version)
import casadi.*
if CL_version == "analytical"
    alpha_mesh = deg2rad([-180, -170, -160, -150, -140, -130, -120, -110, -100, -90, ...
                  -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, ...
                  20, 30, 40, 50, 60, 70, 80, 90, 100, 110, ...
                  120, 130, 140, 150, 160, 170, 180]);

    CL_tail_mesh = [9.82e-02, 4.25e-01, 3.91e-01, 4.27e-01, 4.17e-01, 3.68e-01, 2.80e-01, ...
                    1.58e-01, 1.05e-02, -1.48e-01, -2.98e-01, -4.24e-01, -5.09e-01, ...
                    -5.32e-01, -4.76e-01, -3.29e-01, -1.50e-01, 5.11e-02, -2.45e-01, ...
                    2.14e-01, 5.07e-01, 6.51e-01, 7.03e-01, 6.93e-01, 5.81e-01, 4.44e-01, ...
                    2.71e-01, 8.14e-02, -1.06e-01, -2.74e-01, -4.09e-01, -5.05e-01, ...
                    -5.43e-01, -5.04e-01, -4.72e-01, -4.59e-01, 9.82e-02];

    CL_tail_interp = interpolant('CL_tail_interp', 'bspline', {alpha_mesh}, CL_tail_mesh);

    CL_tail_lookup = CL_tail_interp;

elseif CL_version == "dummy"

    alpha = SX.sym('alpha');
    CL0 = 0.20; CL_alpha = 5.5;
    CL_expr = CL0 + CL_alpha*alpha;
    CL_tail_lookup = Function('CL_tail_lookup', {alpha}, {CL_expr});

else
    warning('CL version unknown')
end
end