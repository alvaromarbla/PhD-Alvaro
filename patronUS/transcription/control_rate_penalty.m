function R = control_rate_penalty(vars, weight)
%CONTROL_RATE_PENALTY Regulariser on the rate of change of the controls.
%   R = weight * sum_i  integral_0^1 (du_i/dtau)^2 dtau,  with u_i the
%   nondimensional control (u/u_scale) and tau = t/tf, approximated by
%   finite differences along the control sequence: u_0, (u_m,0), u_1, ...
%   (midpoint controls included when the method has them). Being an
%   integral, its size does not depend on N or on the method.
%   It removes the flat directions of redundant controls (e.g. tilt vs rps
%   trading the pitching moment) and chattering, at the cost of slightly
%   smoothing bang-bang arcs. weight = 0 disables it.

if weight == 0
    R = 0;
    return
end

[seq, tau] = control_sequence(vars);

dtau = tau(2) - tau(1);
dU   = seq(:,2:end) - seq(:,1:end-1);

R = weight*sumsqr(dU)/dtau;

end
