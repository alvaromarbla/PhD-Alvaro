function [xdot, dbg] = vtolto_dynamics_3GDL(Z, models, params)
% Continuous 3-DOF longitudinal dynamics of the 2Wings aircraft.
% States:   x, h, vx, vh (inertial, h up), theta, q, E (consumed energy)
%           [+ SOC, state of charge, with the Shepherd battery]
% Controls: epsilon1, epsilon2, n1, n2
% xdot has one field per state name.

F = vtolto_forces_3GDL(Z, models, params);

P_bus = (F.P1 + F.P2)/params.prop.eff;          % Electrical power drawn from the battery

xdot.x     = Z.vx;
xdot.h     = Z.vh;
xdot.vx    = F.Fx/params.mass;
xdot.vh    = F.Fh/params.mass;
xdot.theta = Z.q;
xdot.q     = F.M/params.Iyy;
xdot.E     = P_bus;

if strcmp(params.batt.model, 'Shepherd')
    batt = battery_shepherd(P_bus, Z.SOC, params.batt);
    xdot.SOC = -batt.I/params.batt.Q;
    F.I_batt = batt.I;
    F.U_batt = batt.U;
end

F.P_bus = P_bus;
dbg = F;

end
