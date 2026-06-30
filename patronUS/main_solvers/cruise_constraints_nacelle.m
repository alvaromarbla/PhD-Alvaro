function [c, ceq] = cruise_constraints_nacelle(x, models, params, bounds)

% Unpack variables to match the 5-element vector: [V, gamma, alpha, epsilon, n]
V       = x(1);
% gamma = x(2); % Steady cruise implies gamma = 0, which is handled via bounds or optimization rules
alpha   = x(3);
epsilon = x(4);
n       = x(5);

%% Call B-spline models on aero and propulsive forces:
J = V./(n*params.prop.diameter);
phi = alpha + epsilon;

CT = models.CT_lookup(J,phi);
CP = models.CP_lookup(J,phi);
CL = models.CL_lookup(alpha);
CD = models.CD_lookup(alpha);

CD_fuselage = models.CD_fuselage_lookup(alpha);
%% Force calculations

q = 0.5 * params.rho * (V^2) * params.wing_area; % Dynamic pressure

L = CL*q;
D = CD*q;

T =  params.prop.num_engines * CT * params.rho * (n^2) * (params.prop.diameter^4);
P =  params.prop.num_engines * CP * params.rho * (n^3) * (params.prop.diameter^5);


kappa = T/(0.5 * params.rho * params.prop.diameter^2*V^2);
%lambda^4 + 2* sin(phi) *lambda^3 + labda^2 -kappa^2 = 0

coeffs = [1, 2 * sin_phi, 1, 0, -kappa^2];
all_roots = roots(coeffs);

% Filter for the real, positive physical root corresponding to induced velocity
lambda_i = real(all_roots(imag(all_roots) == 0 & real(all_roots) > 0));
alpha_nac = atan(sin(phi)/(lambda_i + cos(phi)));
V_nac2 = V^2 + (lambda_i*V)^2 + 2*cos(phi)*lambda_i*V^2;
q_nac = 0.5 * params.rho * (V_nac2) * params.nacelle_area;

D_nac = CD_fuselage*q_nac;
%% Constraint satisfaction

% Equalities

ceq = [T * cos(phi) - D - D_nac*cos(alpha_nac);  ...                    % Axial force balance
       L - D_nac*sin(alpha_nac )+ T * sin(phi) - (params.mass * params.g)...  % Normal force balance
       ];

% Inequalities

c = [-V/(n*params.prop.diameter)  ; ...
    -CT;...
    -CP; ...
    T - params.prop.T_max_eng*params.prop.num_engines ;...
    P - params.prop.P_max_eng*params.prop.num_engines;
    J-1.2];

end