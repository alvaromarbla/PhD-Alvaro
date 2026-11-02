function defects = collocation_hermite_simpson(F, vars, Fk, h)
%COLLOCATION_HERMITE_SIMPSON Compressed Hermite-Simpson collocation defects (4th order).
%   x_m = (x_k + x_{k+1})/2 + h/8*(f_k - f_{k+1})            (Hermite interpolant)
%   x_{k+1} - x_k - h/6*(f_k + 4*f(x_m, u_m) + f_{k+1}) = 0   (Simpson quadrature)
%   Midpoint controls u_m are decision variables (vars.Um, method.midpoints = true).
%   Fk: state derivatives at the nodes, nx x (N+1)
%   Returns an nx x N matrix of (unscaled) defects.

fk  = Fk(:,1:end-1);
fk1 = Fk(:,2:end);
Xk  = vars.X(:,1:end-1);
Xk1 = vars.X(:,2:end);

Xm = 0.5*(Xk + Xk1) + h/8*(fk - fk1);

F_mid = F.map(vars.N);
fm = F_mid(Xm, vars.Um); % First output: xdot

defects = Xk1 - Xk - h/6*(fk + 4*fm + fk1);

end
