function [P_data] = P_GenerateData_V1(D_data)


% This function generates all data related to propulsion
n_eng = D_data.AC_CONFIGURATION.n_eng;

eta_prop = 0.7;
eta_gen = 0.85; % 15% Related to energy generation
eta_mp = eta_prop*eta_gen;

type_battery = D_data.AC_CONFIGURATION.type_battery;
tau_reserve = 0.10; 
SE_LiFePO4 = 160; %[Wh/kg ] Specific Energy in Joules
SP_motor = 5.82; %[ kW/kg]
SP_ESC = 20; 
m_prop_din = 0.00385;


%% ACTUALIZAR DIAMETRO PROPELLER por el usuario


% 
% Weight_tier.M_ENERGY.SE_LiFePO4 = 160; % Wh/kg Specific Energy  90–160 Wh/kg (320–580 J/g or kJ/kg)
% Weight_tier.M_ENERGY.SE_LiPo = 200; % Wh/kg Specific Energy
% Weight_tier.M_ENERGY.SE_FuelCells = 333; % Wh/kg Specific Energy for fuel cells
% 
% 
% Weight_tier.M_ENERGY.SP_motor = 5.82; % Motor Specific Power kW/kg - TIER 1 
% Weight_tier.M_ENERGY.SP_ESC = 20; % Motor Specific Power kW/kg - TIER 1
% Weight_tier.M_ENERGY.m_prop_din = 0.00385; % Motor Specific Power kW/kg - TIER 1






%% Generate Output 

P_data = load('Prop_data_SANAID');

P_data.type_battery = type_battery;
P_data.tau_reserve = tau_reserve;
P_data.SE_battery = SE_LiFePO4; 
P_data.SP_motor = SP_motor;
P_data.SP_ESC = SP_ESC;
P_data.m_prop_din = m_prop_din;

P_data.eta_prop = eta_prop;
P_data.eta_gen = eta_gen; 
P_data.eta_mp = eta_mp;

P_data.n_eng = n_eng; 

end