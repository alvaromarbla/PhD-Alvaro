function [xdot, dbg] = vtolto_dynamics_3GDL(Z, models, params)
% Continuous 3-DOF longitudinal dynamics of the 2Wings aircraft.
% States:   x, h, vx, vh (inertial, h up), theta, q, E (consumed energy)
% Controls: epsilon1, epsilon2, n1, n2
% xdot has one field per state name.

F = vtolto_forces_3GDL(Z, models, params);

xdot.x     = Z.vx;
xdot.h     = Z.vh;
xdot.vx    = F.Fx/params.mass;
xdot.vh    = F.Fh/params.mass;
xdot.theta = Z.q;
xdot.q     = F.M/params.Iyy;
xdot.E     = (F.P1 + F.P2)/params.prop.eff;

dbg = F;

end
