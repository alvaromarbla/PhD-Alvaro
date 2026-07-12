function [c, ceq] = cruise_constraints_3GDL(x, models, params, bounds)

% Unpack variables to match the 5-element vector: [V, gamma, alpha, epsilon, n]
V       = x(1);
% gamma = x(2); % Steady cruise implies gamma = 0, which is handled via bounds or optimization rules
alpha   = x(3);
epsilon1 = x(4);
n1       = x(5);
epsilon2 = x(6);
n2       = x(7);
delta_e  = x(8);

%% Call B-spline models on aero and propulsive forces:
J = V./(n1*params.prop.diameter);
phi = alpha + epsilon1;

CT = models.CT_lookup(J,phi);
CP = models.CP_lookup(J,phi);
CL = models.CL_lookup(alpha);
CD = models.CD_lookup(alpha);

%% Force calculations

q = 0.5 * params.rho * (V^2) * params.wing_area; % Dynamic pressure

L = CL*q;
D = CD*q;

T =  params.prop.num_engines * CT * params.rho * (n1^2) * (params.prop.diameter^4);
P =  params.prop.num_engines * CP * params.rho * (n1^3) * (params.prop.diameter^5);

%% Constraint satisfaction

% Equalities

ceq = [T * cos(phi) - D;  ...                                  % Axial force balance
       L + T * sin(phi) - (params.mass * params.g)...          % Normal force balance
       ];

% Inequalities

c = [-V/(n1*params.prop.diameter)  ; ...
    -CT;...
    -CP; ...
    T - params.prop.T_max_eng*params.prop.num_engines ;...
    P - params.prop.P_max_eng*params.prop.num_engines;
    J-1.2];

end