function cell = battery_cell(name)
% Cell parameters for the Shepherd model (battery_shepherd, battery_pack).
%   V_oc [V] constant voltage, K [Ohm] polarisation resistance at full
%   charge, R [Ohm] internal resistance, A [V] / B [-] exponential zone,
%   Q [Ah] capacity, mass [kg], V_cut [V] cut-off, I_cont [A] max
%   continuous discharge current.

switch name
    case 'V2'
        % ShepherdModel_V2 (Modules/Battery Models): 2.6 Ah energy cell
        cell = struct('name', name, 'V_oc', 3.7, 'K', 8e-3, 'R', 60e-3, 'A', 0.45, 'B', 4, ...
                      'Q', 2.6, 'mass', 0.045, 'V_cut', 3.0, 'I_cont', NaN);
    case 'P42A'
        % High-power 21700 (Molicel P42A class): R, Q, mass and current
        % from typical datasheet values; the Shepherd shape (V_oc, A, B) is
        % kept from the V2 cell and K is scaled by capacity (2.6/4.2).
        % TYPICAL VALUES, to replace with a fit to the chosen cell
        cell = struct('name', name, 'V_oc', 3.7, 'K', 8e-3*2.6/4.2, 'R', 15e-3, 'A', 0.45, 'B', 4, ...
                      'Q', 4.2, 'mass', 0.070, 'V_cut', 3.0, 'I_cont', 45);
    otherwise
        error('battery_cell:unknown', 'Unknown cell "%s" (V2, P42A).', name);
end

end
