function y = flat_output_spline(t, tf, bc0, bcf)
%FLAT_OUTPUT_SPLINE Quintic polynomial for a flat output between two boundary states.
%   Minimum-order polynomial meeting value, first and second derivative at
%   both ends (6 conditions):
%     bc0 = [y(0),  dy/dt(0),  d2y/dt2(0)]
%     bcf = [y(tf), dy/dt(tf), d2y/dt2(tf)]
%   Returns y.p, y.d, y.dd: the output and its time derivatives at t.

tau = t/tf;

% Boundary conditions in tau = t/tf (d/dtau = tf*d/dt)
b = [bc0(1); bc0(2)*tf; bc0(3)*tf^2; ...
     bcf(1); bcf(2)*tf; bcf(3)*tf^2];

A = [1 0 0 0 0  0;
     0 1 0 0 0  0;
     0 0 2 0 0  0;
     1 1 1 1 1  1;
     0 1 2 3 4  5;
     0 0 2 6 12 20];

coeffs = flipud(A\b)'; % Descending powers for polyval

y.p  = polyval(coeffs, tau);
y.d  = polyval(polyder(coeffs), tau)/tf;
y.dd = polyval(polyder(polyder(coeffs)), tau)/tf^2;

end
