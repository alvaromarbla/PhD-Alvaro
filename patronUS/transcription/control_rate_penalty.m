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

if size(vars.Umn, 2) > 0
    % Interleave nodes and midpoints: u_0, u_m0, u_1, u_m1, ..., u_N
    seq = [reshape([vars.Un(:,1:end-1); vars.Umn], vars.nu, 2*vars.N), vars.Un(:,end)];
else
    seq = vars.Un;
end

dtau = 1/(size(seq, 2) - 1);
dU   = seq(:,2:end) - seq(:,1:end-1);

R = weight*sumsqr(dU)/dtau;

end
