function [D_data,Geo_data] = D_GenerateData_V1
%% Aircraft Configuration

%% PONERLE SCALING FACTOR
D_data = load('AC_CONFIGURATION_SANAID');
n_eng = 4;
D_data.n_eng = n_eng;


%% Geometry Configuration

Geo_data = load('Geo_tier');
Geo_data.Surf_TOT = 1.252259847870206;  % From Body_Geo
Geo_data.f_f = 1.5; 
end
