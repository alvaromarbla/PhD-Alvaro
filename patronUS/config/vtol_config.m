function [params, bounds, models] = vtol_config(versions)
% VTOL take-off segment configuration.
% Aircraft data (mass, geometry, propellers, atmosphere) and models are taken
% from cruise_config so both segments share one aircraft definition; only
% the segment parameters and the trajectory bounds are defined here.

import casadi.*

if versions.AC_version == "2Wings"

    [params, ~, models] = cruise_config(versions);

    % Phase-specific parameters
    params.phase = 'vtolto';
    params.h_to  = 30;                             % [m] Take-off target altitude
    params.x_TO  = 2;                            % [m] Max horizontal displacement at the end of the take-off
    params.V_eps = 0.5;                            % [m/s] Airspeed regularisation at hover

    % Transcription
    params.ocp.N      = 40;                        % Number of collocation intervals
    params.ocp.tf_ref = 10;                        % [s] Reference time (min-time objective scaling)

    % Flat initial guess: built for each tf in the sweep, the feasible one
    % with the lowest objective feeds the first stage
    params.ocp.tf_sweep      = [5 10 15];          % [s]
    params.ocp.guess_dyn_tol = 0.05;               % [-] Max relative dynamics residual of a feasible guess

    tf_init = 15;                                  % [s]

    % States: x0/xf are the boundary conditions (NaN = free)
    % Optional xf_min/xf_max: inequality on the final state
    % init: scalar (constant guess) or [start end] (linear guess)
    % x, vx >= 0: no backwards flight (alpha in vtolto_forces_3GDL assumes u_b >= 0)
    bounds.states.x     = struct('min', 0,            'max', 200,             'x0', 0, 'xf', NaN,         'init', 0, 'xf_max', params.x_TO);
    bounds.states.h     = struct('min', 0,            'max', 1.2*params.h_to,  'x0', 0, 'xf', params.h_to, 'init', [0 params.h_to]);
    bounds.states.vx    = struct('min', 0,         'max', 30,               'x0', 0, 'xf', 0,           'init', 0);
    bounds.states.vh    = struct('min', -5,           'max', 10,               'x0', 0, 'xf', 0,           'init', params.h_to/tf_init);
    bounds.states.theta = struct('min', deg2rad(-20), 'max', deg2rad(20),      'x0', 0, 'xf', 0,           'init', 0);
    bounds.states.q     = struct('min', -0.6,         'max', 0.6,              'x0', 0, 'xf', 0,           'init', 0); % Rad/s
    bounds.states.E     = struct('min', 0,            'max', params.E_batt,    'x0', 0, 'xf', NaN,         'init', [0 0.01*params.E_batt]); % Consumed energy

    % Controls
    bounds.controls.epsilon1 = struct('min', deg2rad(0), 'max', deg2rad(110),      'init', deg2rad(90));
    bounds.controls.epsilon2 = struct('min', deg2rad(0), 'max', deg2rad(110),      'init', deg2rad(90));
    bounds.controls.n1       = struct('min', 10,         'max', params.prop.max_rps, 'init', 50);
    bounds.controls.n2       = struct('min', 10,         'max', params.prop.max_rps, 'init', 50);
    % bounds.controls.deltae = struct('min', deg2rad(-20), 'max', deg2rad(20), 'init', 0); % Not in the tail model yet

    % Final time
    bounds.time.tf = struct('min', 2, 'max', 60, 'init', tf_init);

else
    warning('AC version unknown')

end
