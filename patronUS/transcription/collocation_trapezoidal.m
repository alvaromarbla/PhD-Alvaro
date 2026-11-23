function [defects, Xm] = collocation_trapezoidal(F, vars, Fk, h)
%COLLOCATION_TRAPEZOIDAL Trapezoidal collocation defects (2nd order).
%   x_{k+1} - x_k - h/2*(f_k + f_{k+1}) = 0,  k = 1..N
%   F:  continuous model Function (unused, kept for a common signature)
%   Fk: state derivatives at the nodes, nx x (N+1)
%   Returns an nx x N matrix of (unscaled) defects (and Xm = [], no midpoints).

defects = vars.X(:,2:end) - vars.X(:,1:end-1) - h/2*(Fk(:,1:end-1) + Fk(:,2:end));
Xm = [];   % No midpoints in this scheme

end
