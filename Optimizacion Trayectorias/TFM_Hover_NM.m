%% Code for solving Hovering

clear all
close all



%% Loads and weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;
rho   = 1.2133;

%% Electric
tau = 0.2; 
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency
mbatt = 3; %[kg]
e0    = 720e3; %[J/kg]
E     = e0*mbatt;  %[J] Total energy of the battery packs
N_eng = 2; % Number of engines
D     = 0.8124; % Propeller Diameter [m]


RPMMAX_APC = 150000; % Max RPM of AXI motors 
D_inches = D*1000/25.4; % Diameter in inches
rps_max_APC = RPMMAX_APC/(D_inches); % Max rps of the AXI motors w.r.t. the diameter
rps_max_store = 4687.5;
rps_max = min(rps_max_APC,rps_max_store);

Tmax_Eng = 238.383; % Max Thrust per engine [kN]

Tmax = N_eng*Tmax_Eng;



n_sol = zeros(3,1);
P_sol = zeros(3,1);
t_sol = zeros(3,1);
T_sol = zeros(3,1);
for model = 0:2
   CT = CT_fun(model);
   CP = CP_fun (model);
    
   T = @(V,n) CT(V,n)*N_eng*rho*n^2*D^4;
   P = @(V,n) CP(V,n)*N_eng*rho*n^3*D^5;
    
    %% State function to be solved
Fun = @(n) T(0,n) - W ; % T=W3

n0 = 60;

n_sol(model+1) = fzero(Fun,n0)

T_sol(model+1) = T(0,n_sol(model+1))

P_sol(model+1) = P(0,n_sol(model+1))

t_sol(model+1) = E*eta_m./P_sol(model+1)
end

Tlin = 50;
t_end = 30;

t_vec = linspace(0,t_end,Tlin);
x_vec = zeros(1,Tlin);
h_vec = zeros(1,Tlin);
V_vec = zeros(1,Tlin);
n_vec = n_sol(3,1)*ones(1,Tlin);
gamma_vec = NaN(1,Tlin);
alpha_vec = NaN(1,Tlin);
epsilon_vec = NaN(1,Tlin);
phi_vec = NaN(1,Tlin);
L_vec = NaN(1,Tlin);
T_vec = T_sol(2)*ones(1,Tlin);
D_vec = NaN(1,Tlin);
P_vec = P_sol(2)*ones(1,Tlin);

E_vec = zeros(1,Tlin);
E_vec(1) = E;

for kk = 2:length(t_vec)
E_vec(kk) = E - trapz(P_vec(1:kk)/eta_m)*t_vec(kk)/Tlin;
end

E_vec_perc = (E_vec)/E*100;
% Mission_1 = [t_vec; x_vec ; h_vec; V_vec; n_vec; gamma_vec; alpha_vec ; epsilon_vec; phi_vec; L_vec; T_vec; D_vec; P_vec; E_vec; E_vec_perc]';
% save Mission_1.mat Mission_1
% 
% Mission_3 = [t_vec; x_vec ; h_vec; V_vec; n_vec;gamma_vec; alpha_vec ; epsilon_vec; phi_vec; L_vec; T_vec; D_vec; P_vec; E_vec; E_vec_perc]';
% save Mission_3.mat Mission_3

Mission_10 = [t_vec; x_vec ; h_vec; V_vec; n_vec;gamma_vec; alpha_vec ; epsilon_vec; phi_vec; L_vec; T_vec; D_vec; P_vec; E_vec; E_vec_perc]';
save Mission_10.mat Mission_10

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Propulsive Model

function CT = CT_fun (model)
D     = 0.8124; % Propeller Diameter [m]

if model ==0
% Old model

CT0 =  0.089050757500461;
CT1 = -0.026659496766654;
CT2 = -0.162620845549081;


CT = @(V,n) CT2*V^2/(n^2*D^2)+CT1*V/(n*D)+CT0;

elseif model ==1
    
% Full model

J =@(V,n) V./(n*D);

CT_coeff = [0.0735880531010883, -0.0311758018412727, -0.249744726429543, 0.143084420694372, 0.0261032283758581;
    0, -0.0982459868751664, 0.20127470719351, -0.173738749783189, 0;
    0, 0.156239779501715, 0.0368592048084175, 0 ,0 ;
    0, -0.0478709034281346, 0, 0, 0;
    0, 0, 0, 0, 0];

CT = @(V,n) [1, 0, 0, 0, 0]  * CT_coeff * [1, J(V,n), J(V,n)^2, J(V,n)^3, J(V,n)^4]';

elseif model ==2
    
% App model

CT0 =  0.0735880531010883;
CT1 = -0.0311758018412727;
CT2 = -0.249744726429543;  


CT = @(V,n) CT2*V^2/(n^2*D^2)+CT1*V/(n*D)+CT0;

end

end

%% Propulsive Model

function CP = CP_fun (model)
D     = 0.8124; % Propeller Diameter [m]

if model ==0
% Old model


CP0 = 0.034214580122684;
CP1 = -0.002523708791417;
CP2 = 0.116121898742278;
CP3 = -0.248807360672063;


CP = @(V,n) CP3*V^3/(n^3*D^3)+CP2*V^2/(n^2*D^2)+CP1*V/(n*D)+CP0;

elseif model ==1
    
% Full model

J =@(V,n) V./(n*D);

CP_coeff = [0.0261518307541734, 0.0473735972985378, -0.16267474946046, 0.0247028469343899, 0.0306053713439883;
    0, -0.0762350484603968, 0.148580471912353, -0.0726017200715775, 0;
    0, 0.0897273366920878, 0.0122602815262456, 0, 0;
    0, -0.0486029866039398, 0, 0, 0;
    0, 0, 0, 0, 0];

CP = @(V,n) [1, 0, 0, 0, 0]  * CP_coeff * [1, J(V,n), J(V,n)^2, J(V,n)^3, J(V,n)^4]';

elseif model ==2
    
% App model


CP0 = 0.0261518307541734;
CP1 = 0.0473735972985378;
CP2 = -0.16267474946046;
CP3 = 0.0247028469343899; 


CP = @(V,n) CP3*V^3/(n^3*D^3)+CP2*V^2/(n^2*D^2)+CP1*V/(n*D)+CP0;

end

end


