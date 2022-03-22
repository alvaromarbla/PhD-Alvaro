function [rho_f,rho_fairing,rho_w,rho_HTP,rho_VTP,rho_tb]=CFIIIdensities

% Datos Céfiro III
S_w = 2.200000000000000;
S_w_e = 2.050131939508261;
b_w = 5.225930044268314;
b_w_e = 4.869930044268314;
AR_w = 12.413793103448285;
AR_w_e = 11.568142605375774;
cR_w = 0.420977698010503;
cT_w = 0.420977698010503;

cR_h = 0.259554525531993;
cT_h = 0.259554525531993;
b_h = 1.306482511067079;
S_h = 0.339103448275862;
AR_h = 5.033557046979866;

cR_v = 0.435494170355693;
cT_v = 0.261296502213416;
b_v = 0.478172599050551;
S_v = 0.333186206896552;
AR_v = 1.372500000000000;

% Length Tail boom
l_tb = 2.280000000000000;

% Fuselaje
Surf_TOT = 2.700430097266073; % m^2
Vol_TOT = 0.244461563000000; % m^3
length_fus = 2.552000000000000; % m
w_Area_b_max = 0.356000000000000;
h_Area_b_max = 0.372000000000000;
l_cajon = 1.750000000000000;
Vol_cajon = l_cajon*w_Area_b_max*h_Area_b_max;

% Pesos
W_S_w_e = 9.752;
W_fus = 5.206;
W_S_hv = 3.458;
W_tb = 2.884;
W_fairing = 4.149;

% Recalculo de las nuevas densidades
rho_f = W_fus/Surf_TOT;
rho_fairing = W_fairing/Surf_TOT;
rho_fus = (rho_f + rho_fairing)/2;
rho_w = W_S_w_e/S_w_e;
rho_HTP = W_S_hv/(S_h+S_v);
rho_VTP = W_S_hv/(S_h+S_v);
rho_tb = W_tb/(2*l_tb);
end
