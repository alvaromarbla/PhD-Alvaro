function fh = flight_history(res, params, segment)
% Flight history of a trajectory segment in one standard layout, so that
% the generic post-processing (plot_trajectory_results,
% plot_battery_results, animate_trajectory) does not depend on which
% variables are states, controls or debug outputs of each segment.
%   res:     unpack_ocp_solution output of any segment (VTOL take-off,
%            time-dependent cruise, ...)
%   segment: name used in titles and file names (e.g. 'vtolto', 'cruise')
% Every field is a 1 x n row on the nodes, so the histories of
% consecutive phases can be joined by concatenating their rows (with t
% shifted by the start time of each phase) for a multiphase trajectory.
%
% Each field is taken from the states/controls of res, then from res.dbg,
% and otherwise derived:
%   kinematics  t, x, h, vx, vh, V, gamma, theta, alpha, q
%               (no altitude state: h = params.delta0, steady cruise altitude;
%               no vx/vh: V*cos(gamma), V*sin(gamma))
%   actuators   epsilon1, epsilon2, n1, n2, deltae (NaN if absent)
%   forces      T1, T2, H1, H2, D_nac1, D_nac2, alpha_nac1, alpha_nac2,
%               L_wing, D_wing, L_fus, D_fus, L_tail, D_tail (wind axes),
%               L_aero, D_aero (sum of wing, fuselage and tail)
%   moments     MA_wing, MA_fus, MA_tail, MA, M_Wing, M_Fus, M_Tail,
%               M_Eng1, M_Eng2, M_Nac1, M_Nac2, M (net, = Iyy*dq/dt;
%               ~0 in a trimmed segment)
%   power       P1, P2, P_bus, E, and with the Shepherd battery SOC,
%               I_batt, U_batt
%   meta        segment, tf, n (number of nodes)

d = res.dbg;
n = numel(res.t);
pick = @(varargin) first_field(res, d, n, varargin{:});

fh.segment = segment;
fh.tf = res.tf;
fh.n  = n;
fh.t  = res.t;

%% Kinematics

fh.x     = pick('x');
fh.V     = pick('V');
fh.gamma = pick('gamma');
fh.theta = pick('theta');
fh.alpha = pick('alpha');
if all(isnan(fh.alpha)), fh.alpha = fh.theta - fh.gamma; end
fh.q     = pick('q');
fh.h     = pick('h');
if all(isnan(fh.h))
    h0 = 0;
    if isfield(params, 'delta0'), h0 = params.delta0; end
    fh.h = h0*ones(1, n);
end
fh.vx = pick('vx');
if all(isnan(fh.vx)), fh.vx = fh.V.*cos(fh.gamma); end
fh.vh = pick('vh');
if all(isnan(fh.vh)), fh.vh = fh.V.*sin(fh.gamma); end

%% Actuators

for f = {'epsilon1', 'epsilon2', 'n1', 'n2', 'deltae'}
    fh.(f{1}) = pick(f{1});
end

%% Forces (L_tail/D_tail in wind axes: the cruise names them *_tailwing_wind)

for f = {'T1', 'T2', 'H1', 'H2', 'D_nac1', 'D_nac2', 'alpha_nac1', 'alpha_nac2', 'L_wing', 'D_wing'}
    fh.(f{1}) = pick(f{1});
end
fh.L_fus  = pick('L_fuselage', 'L_fus');
fh.D_fus  = pick('D_fuselage', 'D_fus');
fh.L_tail = pick('L_tailwing_wind', 'L_tail');
fh.D_tail = pick('D_tailwing_wind', 'D_tail');
fh.L_aero = fh.L_wing + fh.L_fus + fh.L_tail;
fh.D_aero = fh.D_wing + fh.D_fus + fh.D_tail;

%% Moments

fh.MA_wing = pick('MA_wing');
fh.MA_fus  = pick('MA_fus', 'MA_fuselage');
fh.MA_tail = pick('MA_tail');
fh.MA      = fh.MA_wing + fh.MA_fus + fh.MA_tail;
for f = {'M_Wing', 'M_Fus', 'M_Tail', 'M_Eng1', 'M_Eng2', 'M_Nac1', 'M_Nac2'}
    fh.(f{1}) = pick(f{1});
end
fh.M = pick('M');
if all(isnan(fh.M))
    fh.M = fh.MA + fh.M_Wing + fh.M_Fus + fh.M_Tail + fh.M_Eng1 + fh.M_Eng2 + fh.M_Nac1 + fh.M_Nac2;
end

%% Power and battery

fh.P1 = pick('P1');
fh.P2 = pick('P2');
fh.P_bus = pick('P_bus');
if all(isnan(fh.P_bus)), fh.P_bus = (fh.P1 + fh.P2)/params.prop.eff; end
fh.E = pick('E');
for f = {'SOC', 'I_batt', 'U_batt'}
    fh.(f{1}) = pick(f{1});
end

end

function v = first_field(res, d, n, varargin)
% First of the names found in res (states/controls) or res.dbg, as a row;
% NaN row if none exists
for k = 1:numel(varargin)
    name = varargin{k};
    if isfield(res, name) && isnumeric(res.(name)) && numel(res.(name)) == n
        v = reshape(res.(name), 1, n);
        return
    end
    if isfield(d, name)
        v = reshape(d.(name), 1, n);
        return
    end
end
v = NaN(1, n);
end
