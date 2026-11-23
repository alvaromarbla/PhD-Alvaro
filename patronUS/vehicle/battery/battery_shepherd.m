function out = battery_shepherd(P_bus, SOC, batt)
% Shepherd battery model at pack level (ShepherdModel_V2 / report form),
% valid for any power history P_bus(t):
%     U = E0 - R*I - K*Q/(Q - q)*I + A*exp(-B*q/Q),     q = (1 - SOC)*Q
%       = E0 - (R + K/SOC)*I + A*exp(-B*(1 - SOC))
% The power balance U*I = P_bus is a quadratic in the current,
%     (R + K/SOC)*I^2 - U_oc*I + P_bus = 0,   U_oc = E0 + A*exp(-B*(1 - SOC)),
% whose smaller root is the physical current. It is written in the
% cancellation-free form I = 2*P_bus/(U_oc + sqrt(disc)).
% disc = U_oc^2 - 4*(R + K/SOC)*P_bus >= 0 means the battery can deliver
% P_bus; inside the expression it is replaced by a smooth max(disc, 0)
% (width batt.disc_eps [V^2]) so that the current stays real during the
% NLP iterations. The constraint disc >= 0 is imposed separately.
% Works with numbers and with CasADi symbols.
%   P_bus [W], SOC [-], batt: pack parameters (vtol_config)
%   out.I [A], out.U [V], out.U_oc [V], out.disc [V^2]

a      = batt.R + batt.K./SOC;                  % Internal + polarisation resistance [Ohm]
U_oc   = batt.E0 + batt.A.*exp(-batt.B.*(1 - SOC));
disc   = U_oc.^2 - 4*a.*P_bus;
disc_s = 0.5*(disc + sqrt(disc.^2 + batt.disc_eps^2));

out.I    = 2*P_bus./(U_oc + sqrt(disc_s));
out.U    = U_oc - a.*out.I;
out.U_oc = U_oc;
out.disc = disc;

end
