function [seq, tau] = control_sequence(vars)
%CONTROL_SEQUENCE Nondimensional controls in time order, with their tau = t/tf.
%   Nodes only, or nodes and midpoints interleaved (u_0, u_m0, u_1, ...,
%   u_N) when the method has midpoint controls. Used by the control-rate
%   regulariser and the actuator rate limits.

if size(vars.Umn, 2) > 0
    seq = [reshape([vars.Un(:,1:end-1); vars.Umn], vars.nu, 2*vars.N), vars.Un(:,end)];
    tau = reshape([vars.tau(1:end-1); vars.tau_mid], 1, []);
    tau = [tau, vars.tau(end)];
else
    seq = vars.Un;
    tau = vars.tau;
end

end
