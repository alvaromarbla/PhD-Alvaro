function guess = cruise_soc_guess(trim, P_bus, params, bounds)
% Initial guess of the time-dependent cruise from a steady cruise point:
% the trim variables are held constant (trim: struct with one value per
% field of bounds.controls, e.g. the steady optimum) and the Shepherd
% battery is discharged at the constant bus power P_bus [W] from SOC0 to
% SOC_f (explicit integration with a 1 s step), which gives tf, the range
% and the energy drawn.

batt  = params.batt;
SOC_f = bounds.states.SOC.xf;
dt    = 1;                                           % [s]
t_max = bounds.time.tf.max;

n   = floor(t_max/dt) + 1;
SOC = zeros(1, n);
SOC(1) = batt.SOC0;
k = 1;
while SOC(k) > SOC_f && k < n
    op = battery_shepherd(P_bus, SOC(k), batt);
    if op.disc < 0
        break                                        % Power not deliverable any more
    end
    SOC(k+1) = SOC(k) - op.I/batt.Q*dt;
    k = k + 1;
end
SOC = SOC(1:k);
t   = (0:k-1)*dt;

% Last step interpolated to land on SOC_f
if SOC(end) < SOC_f && k > 1
    t(end)   = t(end-1) + dt*(SOC(end-1) - SOC_f)/(SOC(end-1) - SOC(end));
    SOC(end) = SOC_f;
end

V = trim.V;
if isfield(trim, 'gamma'), V = V*cos(trim.gamma); end

guess.tf  = t(end);
guess.t   = t;
guess.SOC = SOC;
guess.x   = V*t;
guess.E   = P_bus*t;

names = fieldnames(bounds.controls);
for i = 1:numel(names)
    guess.(names{i}) = trim.(names{i})*ones(size(t));
end

end
