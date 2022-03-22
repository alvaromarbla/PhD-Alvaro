function [R_data] = R_GenerateData_V1(W_data,Geo_data)

global g
% Generation of Performance Data for the SAA

R_data = load('Performance');


SafetyF = 1.2; % Safety factor for stall speed


m_TOW     = W_data.m_TOW; % A/C MTOW
S_w1      = Geo_data.S_w1_ph;% Reference surface
CL_max_w1 = R_data.Performance.CL_max_w1; % Maximum operative Lift coefficient
rho       = R_data.Performance.rho; % Air density at h = h_climb m


V_stall = sqrt(2*m_TOW*g/(rho*S_w1*CL_max_w1)); % Stall speed

Vh_min                     = V_stall*SafetyF; % Minimum operating speed
Vh_max                     = Vh_min*1.2; % Max speed. 1.2 is not a safety factor, can be changed


R_data.Performance.V_stall = V_stall;

R_data.Performance.V_min   = Vh_min;
R_data.Vh_vect             = linspace(Vh_min,Vh_max);


R_data.Vv_min              = 0.2; %[m/s] Minimum speed for the VTOL
R_data.Performance.Climbh  = 50; % [m] Climb height

end