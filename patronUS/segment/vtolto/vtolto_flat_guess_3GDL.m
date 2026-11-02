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

n_pts = numel(t);
[epsilon, n1, n2, P] = deal(zeros(1, n_pts));

for k = 1:n_pts
    Z = struct('x', xs.p(k), 'h', hs.p(k), 'vx', xs.d(k), 'vh', hs.d(k), ...
               'theta', ths.p(k), 'q', ths.d(k), 'E', 0, ...
               'epsilon1', pi/2, 'epsilon2', pi/2, 'n1', mean(C.n1.init), 'n2', mean(C.n2.init));

    % Aero forces and moments only depend on the states
    F = vtolto_forces_3GDL(Z, models, params);
    V     = full(F.V);
    alpha = full(F.alpha);
    L_aero = full(F.L_aero);
    D_aero = full(F.D_aero);
    cg = Z.vx/V;
    sg = Z.vh/V;

    % Thrust vector required by the force balance (inertial axes)
    Tx = params.mass*xs.dd(k) - (-D_aero*cg - L_aero*sg);
    Th = params.mass*(hs.dd(k) + params.g) - (-D_aero*sg + L_aero*cg);
    T_tot = hypot(Tx, Th);

    epsilon(k) = min(max(atan2(Th, Tx) - Z.theta, C.epsilon1.min), C.epsilon1.max);

    % Differential thrust to trim the moment:
    % (T1-T2)*sin(eps)*xtw - (T1+T2)*cos(eps)*ztw = Iyy*thetaddot - M_aero
    M_aero = full(F.MA + F.M_Wing + F.M_Fus + F.M_Tail);
    M_req  = params.Iyy*ths.dd(k) - M_aero;
    dT = 0;
    if abs(sin(epsilon(k))) > 0.2
        dT = (M_req + T_tot*cos(epsilon(k))*params.geo.ztw)/(sin(epsilon(k))*params.geo.xtw);
        dT = min(max(dT, -0.9*T_tot), 0.9*T_tot);
    end

    phi = alpha + epsilon(k);
    n1(k) = thrust_to_rps((T_tot + dT)/2, V, phi, C.n1, models, params);
    n2(k) = thrust_to_rps((T_tot - dT)/2, V, phi, C.n2, models, params);

    % Power with the inverted controls
    Z.epsilon1 = epsilon(k);  Z.epsilon2 = epsilon(k);
    Z.n1 = n1(k);             Z.n2 = n2(k);
    F = vtolto_forces_3GDL(Z, models, params);
    P(k) = full(F.P1 + F.P2);
end

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

function n = thrust_to_rps(T, V, phi, cfg, models, params)
% Solve T = Ne*CT(J(n), phi)*rho*n^2*D^4 for n by fixed-point iteration
D  = params.prop.diameter;
Ne = params.prop.num_engines;

n = mean(cfg.init);
for it = 1:20
    CT = full(models.CT_lookup(V/(n*D), phi));
    if CT <= 0
        n = cfg.max;
        return
    end
    n = min(max(sqrt(max(T, 0)/(Ne*CT*params.rho*D^4)), cfg.min), cfg.max);
end
end
