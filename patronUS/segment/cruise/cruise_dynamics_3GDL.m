function [xdot, dbg] = cruise_dynamics_3GDL(Z, models, params)
% Time evolution of the quasi-steady cruise (2Wings): the trim variables
% (V, theta, epsilon1/2, n1/2, ...) are controls at every node and only the
% range, the bus energy and, with the Shepherd battery, the state of charge
% are integrated:
%     dx/dt   = V*cos(gamma)
%     dE/dt   = P_bus = (P1 + P2)/eta
%     dSOC/dt = -I(P_bus, SOC)/Q          (battery_shepherd)
% P1, P2 from the propeller power model, as in cruise_constraints_3GDL.

V     = Z.V;
gamma = Z.gamma;
alpha = Z.theta - gamma;

D  = params.prop.diameter;
J1 = V/(Z.n1*D);
J2 = V/(Z.n2*D);
CP1 = models.CP_lookup(J1, alpha + Z.epsilon1);
CP2 = models.CP_lookup(J2, alpha + Z.epsilon2);
P1 = params.prop.num_engines*CP1*params.rho*(Z.n1^3)*(D^5);
P2 = params.prop.num_engines*CP2*params.rho*(Z.n2^3)*(D^5);

P_bus = (P1 + P2)/params.prop.eff;

xdot.x = V*cos(gamma);
xdot.E = P_bus;

dbg.P_bus = P_bus;
if strcmp(params.batt.model, 'Shepherd')
    batt = battery_shepherd(P_bus, Z.SOC, params.batt);
    xdot.SOC   = -batt.I/params.batt.Q;
    dbg.I_batt = batt.I;
    dbg.U_batt = batt.U;
end

end
