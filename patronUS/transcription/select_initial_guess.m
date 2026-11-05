function [guess, sweep] = select_initial_guess(guess_fcn, tf_sweep, objective_fcn, dynamics_fcn, constraint_fcn, bounds, models, params)
%SELECT_INITIAL_GUESS Sweep the duration of a flat-output guess and keep the best one.
%   For each tf in tf_sweep builds guess_fcn(bounds, models, params, tf),
%   evaluates the objective on it and checks it with check_ocp_guess.
%   Returns the feasible guess with the lowest objective; if none is
%   feasible, the least infeasible one (violation + dynamics residual),
%   with a warning. guess_fcn returning [] (config guess) returns [].
%   sweep is a struct array with tf, obj and the check of every candidate.

guess = [];
sweep = struct('tf', {}, 'obj', {}, 'feasible', {}, 'bound_viol', {}, ...
               'bc_viol', {}, 'path_viol', {}, 'dyn_res', {}, 'trim_viol', {});
candidates = cell(1, numel(tf_sweep));

for k = 1:numel(tf_sweep)
    g = guess_fcn(bounds, models, params, tf_sweep(k));
    if isempty(g)
        return
    end
    chk = check_ocp_guess(g, dynamics_fcn, constraint_fcn, bounds, models, params);
    candidates{k} = g;
    sweep(k).tf  = tf_sweep(k);
    sweep(k).obj = full(objective_fcn(g, models, params));
    for f = {'feasible', 'bound_viol', 'bc_viol', 'path_viol', 'dyn_res', 'trim_viol'}
        sweep(k).(f{1}) = chk.(f{1});
    end
end

%% Report

fprintf('\nInitial guess sweep:\n');
fprintf('   tf [s]   objective   feasible   bounds     BCs        path       dyn_res    trim\n');
for k = 1:numel(sweep)
    fprintf('  %6.2f   %9.4f   %8d   %8.2e   %8.2e   %8.2e   %8.2e   %8.2e\n', sweep(k).tf, sweep(k).obj, ...
        sweep(k).feasible, sweep(k).bound_viol, sweep(k).bc_viol, sweep(k).path_viol, sweep(k).dyn_res, sweep(k).trim_viol);
end

%% Selection

obj      = [sweep.obj];
feasible = [sweep.feasible];

if any(feasible)
    obj(~feasible) = Inf;
    [~, best] = min(obj);
else
    viol = max([[sweep.bound_viol]; [sweep.bc_viol]; [sweep.path_viol]], [], 1);
    [~, best] = min(viol + [sweep.dyn_res]);
    warning('select_initial_guess:infeasible', ...
        'No feasible guess in the sweep, using the least infeasible one (tf = %g s).', sweep(best).tf);
end

guess = candidates{best};
fprintf('Selected guess: tf = %g s (objective %.4f)\n\n', sweep(best).tf, sweep(best).obj);

end
