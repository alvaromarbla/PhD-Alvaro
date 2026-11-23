function batt = battery_pack(cell, Ns, Np)
% Shepherd pack parameters from a cell (battery_cell) and its Ns x Np layout:
% voltages scale with Ns, resistances with Ns/Np, capacity with Np.
% Q is returned in As, as used by the dynamics (dSOC/dt = -I/Q).

batt.model  = 'Shepherd';
batt.cell   = cell.name;
batt.Ns     = Ns;
batt.Np     = Np;
batt.E0     = cell.V_oc*Ns;
batt.A      = cell.A*Ns;
batt.B      = cell.B;
batt.R      = cell.R*Ns/Np;
batt.K      = cell.K*Ns/Np;
batt.Q      = cell.Q*Np*3600;          % [As]
batt.V_cut  = cell.V_cut*Ns;           % [V]
batt.I_cont = cell.I_cont*Np;          % [A] Max continuous pack current
batt.mass   = cell.mass*Ns*Np;         % [kg]

end
