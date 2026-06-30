function [t, y] = DAESolver(odefun, tspan, y0, options)
    M = diag([1 1 1 0]); % Mass matrix for DAE system
    options = odeset(options, 'Mass', M);
    [t, y] = ode15s(odefun, tspan, y0, options);
end