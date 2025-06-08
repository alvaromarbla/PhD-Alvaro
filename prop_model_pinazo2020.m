function [CT, CP] = prop_model_pinazo2020 (V,n,phi)


CT0 =  0.0735880531010883;
CT1 = -0.0311758018412727;
CT2 = -0.249744726429543;  


CP0 = 0.0261518307541734;
CP1 = 0.0473735972985378;
CP2 = -0.16267474946046;
CP3 = 0.0247028469343899;


%% Power
CP =  CP3*V^3/(n^3*D^3)+CP2*V^2/(n^2*D^2)+CP1*V/(n*D)+CP0;

P =  CP(V,n)*N_eng*rho*n^3*D^5;

%% Thrust
CT = @(V,n) CT2*V^2/(n^2*D^2)+CT1*V/(n*D)+CT0;

T = @(V,n) CT(V,n)*N_eng*rho*n^2*D^4;