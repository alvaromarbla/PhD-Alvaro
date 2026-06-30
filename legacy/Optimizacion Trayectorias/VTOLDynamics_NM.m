function dz = VTOLDynamics_NM(z,u)
% dz = rocketDynamics(z,u)
%
% The basic dynamics and drag coefficient data are from the paper:
%
%   "Drag-law Effects in the Goddard Problem"
%   P. Tsiotras, H. Kelley, H.Kelley    1991
%
% INPUTS:
%   z = [3,n] = [h; v; m] = state vector
%   u = [1,n] = [T] = control = thrust
%

h = z(1,:);   %Height
V = z(2,:);   %Velocity
E = z(3,:);   %Energy

n = u;   %RPS

%%%% Density of air:
rho = 1.2133; % Air density [kg/m^3]
D   = 0.8128; % Propeller Diameter [m]
N_eng = 2; % Number of engines

%% Aerodynamic model
alpha_V = -90*pi/180; % [º] Because of Vertical TO


p_CD_ac = [0.08034; 0.0179749;3.715949;-0.852857;-2.376796;1.049514;0.53177; -0.485949; -0.0586846;0.0778734];

%% Loads and weight
mTOW = 16.6; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
S = 0.5008; % Reference surface [m^2]


%% Electric
tau = 0.2;
eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency

% Full model w.r.t angle of attack

CD = @(alpha) p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + p_CD_ac(4)*alpha.^3 ...
    + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);


%% Propulsive Model

CT0 =  0.0735880531010883;
CT1 = -0.0311758018412727;
CT2 = -0.249744726429543; 

CP0 = 0.0261518307541734;
CP1 = 0.0473735972985378;
CP2 = -0.16267474946046;
CP3 = 0.0247028469343899;

%% Thrust
CT = @(V,n) CT2*V.^2./(n.^2*D^2)+CT1*V./(n*D)+CT0;

T = @(V,n) CT(V,n).*N_eng*rho.*n.^2*D^4;

CP = @(V,n) CP3*V.^3./(n.^3*D^3)+CP2*V.^2./(n.^2*D^2)+CP1.*V./(n*D)+CP0;
P = @(V,n)  CP(V,n).*N_eng*rho.*n.^3*D^5;

%%%% Complete system of equations
dh = V;  %vertical velocity

Jfactor = V./(n*D);

Tfactor = T(V,n);
Dfactor = 1/2*rho*V.^2*S*CD(alpha_V);

dv = T(V,n)./mTOW - 1/2*rho*V.^2*S*CD(alpha_V)/mTOW- g;   %vertical acceleration
dE = P(V,n)/(eta_m);   %mass rate

dz = [dh;dv;dE];

end