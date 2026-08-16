% Generation of AERODYNAMIC LIFT model for Fuselage based on results by
% Samuel Torres (2022) For AC with 2 wings. This only accounts the fuselage
% to add it up to other components

function CL_fuselage_lookup = cl_fuselage(CL_version)
import casadi.*

if CL_version == "analytical"

alpha_mesh = deg2rad([-180, -170, -160, -150, -140, -130, -120, -110, -100, -90, ...
    -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, ...
    20, 30, 40, 50, 60, 70, 80, 90, 100, 110, ...
    120, 130, 140, 150, 160, 170, 180]);

CL_fuselage_2wing_mesh = [-1.28e-02, 2.28e-02, 1.43e-01, 3.16e-01, 4.62e-01, 5.53e-01, 5.66e-01, ...
    4.37e-01, 1.87e-01, 5.04e-02, -7.69e-02, -2.33e-01, -5.56e-01, ...
    -6.32e-01, -6.15e-01, -5.04e-01, -3.28e-01, -1.56e-01, -3.68e-02, ...
    8.80e-03, 8.69e-02, 2.54e-01, 4.40e-01, 6.38e-01, 7.56e-01, 5.27e-01, ...
    2.33e-01, 8.31e-02, -7.31e-02, -2.32e-01, -4.33e-01, -6.55e-01, ...
    -6.70e-01, -5.74e-01, -3.73e-01, -9.69e-02, -1.28e-02];

CL_fuselage_interp = interpolant('CL_fuselage_interp','bspline',{alpha_mesh},CL_fuselage_2wing_mesh);

% Symbolic translation using CASADI for OCP:


CL_fuselage_lookup = CL_fuselage_interp;

else
    warning('CL version unknown')

end

end