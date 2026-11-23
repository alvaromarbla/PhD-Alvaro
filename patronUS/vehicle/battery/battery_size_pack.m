function [Np, info] = battery_size_pack(cell, Ns, P_req, margin)
% Smallest number of parallel strings Np such that, with a full pack
% (Shepherd model, battery_shepherd):
%   1. the power deliverable above the cut-off voltage is at least
%      (1 + margin)*P_req:   V_cut*(U_oc - V_cut)/(R + K) >= (1 + margin)*P_req
%   2. drawing P_req, the current per cell stays within I_cont.
% P_req [W] is the peak bus power the pack must supply.
% info: pack power at cut-off, current and terminal voltage at P_req.

for Np = 1:200
    batt = battery_pack(cell, Ns, Np);
    batt.disc_eps = 0;
    a    = batt.R + batt.K;
    U_oc = batt.E0 + batt.A;
    P_cut = batt.V_cut*(U_oc - batt.V_cut)/a;
    op    = battery_shepherd(P_req, 1, batt);
    ok_power   = P_cut >= (1 + margin)*P_req;
    ok_current = isnan(cell.I_cont) || (op.disc >= 0 && op.I/Np <= cell.I_cont);
    if ok_power && ok_current
        break
    end
end

info.P_cut  = P_cut;
info.I      = op.I;
info.I_cell = op.I/Np;
info.U      = op.U;
info.mass   = batt.mass;
info.limit  = '';
if Np > 1
    prev = battery_pack(cell, Ns, Np - 1); prev.disc_eps = 0;
    P_cut_prev = prev.V_cut*(prev.E0 + prev.A - prev.V_cut)/(prev.R + prev.K);
    if P_cut_prev < (1 + margin)*P_req
        info.limit = 'power at cut-off';
    else
        info.limit = 'cell current';
    end
end

end
