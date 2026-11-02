function slack = ocp_bound_slack(W_hist, vars, bounds)
%OCP_BOUND_SLACK Bound slack per variable name for a transcribed OCP.
%   Same output as computeBoundSlack (so plotConvergenceMetrics can use it),
%   but each variable spans several nodes, so the tightest node is reported.
%   W_hist columns follow vars.w = [X(:); U(:); Um(:); tf].

n_iter = size(W_hist, 2);
nX  = vars.nx*(vars.N+1);
nU  = vars.nu*(vars.N+1);
nUm = vars.nu*size(vars.Um, 2);

X_hist  = reshape(W_hist(1:nX, :),                 vars.nx, vars.N+1, n_iter);
U_hist  = reshape(W_hist(nX+1:nX+nU, :),           vars.nu, vars.N+1, n_iter);
Um_hist = reshape(W_hist(nX+nU+1:nX+nU+nUm, :),    vars.nu, [],       n_iter);
tf_hist = W_hist(end, :);

slack = struct();
for i = 1:vars.nx
    name = vars.state_names{i};
    slack.(name) = node_slack(reshape(X_hist(i,:,:), [], n_iter), bounds.states.(name));
end
for i = 1:vars.nu
    name = vars.control_names{i};
    vals = [reshape(U_hist(i,:,:), [], n_iter); reshape(Um_hist(i,:,:), [], n_iter)];
    slack.(name) = node_slack(vals, bounds.controls.(name));
end
slack.tf = node_slack(tf_hist, bounds.time.tf);

end

function s = node_slack(vals, cfg)
% vals: n_nodes x n_iter
s.lower      = min(vals - cfg.min, [], 1);
s.upper      = min(cfg.max - vals, [], 1);
s.tightest   = min(s.lower, s.upper);
s.normalized = s.tightest ./ max(cfg.max - cfg.min, eps);
end
