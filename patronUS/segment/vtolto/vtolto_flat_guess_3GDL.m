function guess = vtolto_flat_guess_3GDL(bounds, models, params, tf)
% Initial guess for the VTOL take-off by differential flatness.
% Flat outputs x(t), h(t), theta(t) are quintic splines meeting the position,
% velocity and (zero) acceleration boundary conditions of vtol_config, with
% x(tf) = params.x_TO. tf is the guess duration (default bounds.time.tf.init;
% select_initial_guess sweeps it). The other states follow
% by differentiation, and the controls by inverting the force and moment
% balance point by point. Only a guess, so it assumes:
%   - epsilon1 = epsilon2, along the required thrust direction
%   - H = 0 and no nacelle drag
%   - n1, n2 split the thrust so the engines trim the pitching moment
%     (equal split when the tilt is too close to horizontal for that)
% guess has the layout of unpack_ocp_solution (t, tf, <state>, <control>),
% so it is passed to set_ocp_initial_guess.

S = bounds.states;
C = bounds.controls;

if nargin < 4
    tf = bounds.time.tf.init;
end
t = linspace(0, tf, 101);

%% Flat outputs

xs  = flat_output_spline(t, tf, [S.x.x0,     S.vx.x0, 0], [params.x_TO, S.vx.xf, 0]);
hs  = flat_output_spline(t, tf, [S.h.x0,     S.vh.x0, 0], [S.h.xf,      S.vh.xf, 0]);
ths = flat_output_spline(t, tf, [S.theta.x0, S.q.x0,  0], [S.theta.xf,  S.q.xf,  0]);

%% Controls by inversion
% Vectorised over all points: the forces are evaluated through one mapped
% CasADi Function (numeric DM arithmetic point by point is very slow)

n_pts = numel(t);
F_pts = forces_function(models, params).map(n_pts);
ones_pts = ones(1, n_pts);

% Aero forces and moments only depend on the states (dummy controls)
F = eval_forces(F_pts, xs.d, hs.d, ths.p, pi/2*ones_pts, pi/2*ones_pts, ...
                mean(C.n1.init)*ones_pts, mean(C.n2.init)*ones_pts);
cg = xs.d./F.V;
sg = hs.d./F.V;

% Thrust vector required by the force balance (inertial axes)
Tx = params.mass*xs.dd - (-F.D_aero.*cg - F.L_aero.*sg);
Th = params.mass*(hs.dd + params.g) - (-F.D_aero.*sg + F.L_aero.*cg);
T_tot = hypot(Tx, Th);

epsilon = min(max(atan2(Th, Tx) - ths.p, C.epsilon1.min), C.epsilon1.max);

% Thrust split to trim the moment, with the rotor arms of M_Eng1/M_Eng2 in
% vtolto_forces_3GDL (rotor 1 at xw, rotor 2 at -xtw, both at ztw; H = 0):
% T1*sin(eps)*xw - T2*sin(eps)*xtw - (T1+T2)*cos(eps)*ztw = Iyy*thetaddot - M_aero
g = params.geo;
M_aero = F.MA + F.M_Wing + F.M_Fus + F.M_Tail;
M_req  = params.Iyy*ths.dd - M_aero;
T1   = T_tot/2;
trim = abs(sin(epsilon)) > 0.2;
T1(trim) = (M_req(trim) + T_tot(trim).*(cos(epsilon(trim))*g.ztw + sin(epsilon(trim))*g.xtw)) ...
           ./(sin(epsilon(trim))*(g.xw + g.xtw));
T1 = min(max(T1, 0.05*T_tot), 0.95*T_tot);
T2 = T_tot - T1;

phi = F.alpha + epsilon;
n1 = thrust_to_rps(T1, F.V, phi, C.n1, models, params);
n2 = thrust_to_rps(T2, F.V, phi, C.n2, models, params);

% Power with the inverted controls
F = eval_forces(F_pts, xs.d, hs.d, ths.p, epsilon, epsilon, n1, n2);
P = F.P1 + F.P2;

%% Output

guess.t  = t;
guess.tf = tf;

guess.x     = xs.p;
guess.h     = hs.p;
guess.vx    = xs.d;
guess.vh    = hs.d;
guess.theta = ths.p;
guess.q     = ths.d;
guess.E     = cumtrapz(t, P/params.prop.eff);

guess.epsilon1 = epsilon;
guess.epsilon2 = epsilon;
guess.n1       = n1;
guess.n2       = n2;

end

function F = forces_function(models, params)
% vtolto_forces_3GDL as a CasADi Function of [vx; vh; theta; epsilon1; epsilon2; n1; n2],
% returning the fields used here stacked in one column
import casadi.*
z = MX.sym('z', 7);
Z = struct('vx', z(1), 'vh', z(2), 'theta', z(3), ...
           'epsilon1', z(4), 'epsilon2', z(5), 'n1', z(6), 'n2', z(7));
Fs = vtolto_forces_3GDL(Z, models, params);
out = cellfun(@(name) Fs.(name), force_names(), 'UniformOutput', false);
F = Function('F_forces', {z}, {vertcat(out{:})});
end

function names = force_names()
names = {'V', 'alpha', 'L_aero', 'D_aero', 'MA', 'M_Wing', 'M_Fus', 'M_Tail', 'P1', 'P2'};
end

function F = eval_forces(F_pts, vx, vh, theta, epsilon1, epsilon2, n1, n2)
% Evaluate the mapped Function at all points, one row per field
vals  = full(F_pts([vx; vh; theta; epsilon1; epsilon2; n1; n2]));
names = force_names();
for i = 1:numel(names)
    F.(names{i}) = vals(i,:);
end
end

function n = thrust_to_rps(T, V, phi, cfg, models, params)
% Solve T = Ne*CT(J(n), phi)*rho*n^2*D^4 for n by fixed-point iteration (all points at once)
D  = params.prop.diameter;
Ne = params.prop.num_engines;
CT_pts = models.CT_lookup.map(numel(T));

n = mean(cfg.init)*ones(size(T));
no_thrust = false(size(T)); % CT <= 0: the propeller cannot give the thrust, saturate
for it = 1:20
    CT = full(CT_pts(V./(n*D), phi));
    no_thrust = no_thrust | CT <= 0;
    n = min(max(sqrt(max(T, 0)./(Ne*max(CT, eps)*params.rho*D^4)), cfg.min), cfg.max);
    n(no_thrust) = cfg.max;
end
end
