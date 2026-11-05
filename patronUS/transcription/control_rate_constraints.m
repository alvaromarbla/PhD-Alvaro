function [c, rows] = control_rate_constraints(vars, bounds)
%CONTROL_RATE_CONSTRAINTS Actuator rate limits |du/dt| <= rate_max.
%   For every control with bounds.controls.<name>.rate_max, the rate between
%   consecutive control points (nodes, and midpoints when the method has
%   them) is limited, normalised by the limit:
%       (du/dt)/rate_max - 1 <= 0   and   -(du/dt)/rate_max - 1 <= 0
%   c stacks all of them (c <= 0). rows describes each side for
%   report_constraint_activity: name, index into c, the symbolic rate, the
%   limit and the tau of each interval centre.

import casadi.*

[seq, tau] = control_sequence(vars);
seq = diag(vars.u_scale)*seq;          % Physical units
dt  = (tau(2) - tau(1))*vars.tf;
tau_mid = 0.5*(tau(1:end-1) + tau(2:end));

c    = MX(0, 1);
rows = struct('name', {}, 'idx', {}, 'rate', {}, 'limit', {}, 'tau', {});

for i = 1:vars.nu
    name = vars.control_names{i};
    cfg  = bounds.controls.(name);
    if ~isfield(cfg, 'rate_max') || isinf(cfg.rate_max)
        continue
    end
    rate = (seq(i,2:end) - seq(i,1:end-1))/dt;
    n = size(rate, 2);
    for side = [1, -1]
        idx = size(c, 1) + (1:n);
        c   = [c; (side*rate/cfg.rate_max - 1)']; %#ok<AGROW>
        if side > 0
            label = sprintf('d%s/dt <= %g', name, cfg.rate_max);
        else
            label = sprintf('d%s/dt >= %g', name, -cfg.rate_max);
        end
        rows(end+1) = struct('name', label, 'idx', idx, 'rate', rate, ...
                             'limit', side*cfg.rate_max, 'tau', tau_mid); %#ok<AGROW>
    end
end

end
