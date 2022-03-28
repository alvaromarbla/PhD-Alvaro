%% Code for solving Max Range for a given Energy
clear all
close all
mbatt = 7; %[kg]
e0    = 720e3; %[J/kg]

E     = e0*mbatt;  %[J] Total energy of the battery packs

N_eng = 2; % Number of engines
D     = 0.7112; % Propeller Diameter [m]


CT0 =  0.089050757500461;
CT1 = -0.026659496766654;
CT2 = -0.162620845549081; 


CP0 = 0.034214580122684;
CP1 = -0.002523708791417;
CP2 = 0.116121898742278;
CP3 = -0.248807360672063;

mTOW = 21.39; % Maximum T-O Mass [kg]
g    = 9.81; % Gravity [m/s^2]
W     = mTOW*g;

S_ref = 0.430279179101; % Reference surface [m^2]
rho   = 1.225;

tau = 0.2; 
CD0 = 0.035682195723975;
CD2 = 0.054209627025009;

eta_m = 0.88*0.98*0.95; % Engine and electrical efficiency




%% Power
CP = @(V,n) CP3*V^3/(n^3*D^3)+CP2*V^2/(n^2*D^2)+CP1*V/(n*D)+CP0;

P = @(V,n) CP(V,n)*N_eng*rho*n^3*D^5;


%% Range

xf = @(V,n) eta_m*E*V/((1-tau)*P(V,n))*1e-3;

%% Generate numerical Matrix
Nmat = 20;
Vlin = linspace(20,40,Nmat);
nlin = linspace(40,100, Nmat);

% Preallocate matrix for speed
xfmat = zeros(Nmat,Nmat);
% Loop in speed
for  ii = 1: Nmat
    % Loop in revolutions
    for jj = 1:Nmat
        
        xfmat(ii,jj) = xf(Vlin(ii),nlin(jj));
        
    end
end

%% Calculate restriction T = D 

ncon = @(V) -(1/2)*(CT1*N_eng*S_ref*D*V^2*rho-sqrt(-4*CT0*CT2*D^2*N_eng^2*S_ref^2*V^4*rho^2+CT1^2*D^2*N_eng^2*S_ref^2*V^4*rho^2+2*CD0*CT0*N_eng*S_ref^3*V^4*rho^2+...
    8*CT0*N_eng*S_ref*W^2*CD2))/(CT0*N_eng*S_ref*rho*V*D^2);
% Preallocate matrix for speed
nconlin = zeros(Nmat,1);
for kk = 1:Nmat
    nconlin (kk) = ncon(Vlin(kk));

end

%% Eliminate zeros from indeterminations (points where CP = 0)

threshold = 1e2; % Assume Range < 100 km 
xfmat(xfmat > threshold) = 0;
xfmat(xfmat < 0) = 0;


%% Plot contours
N_contour_lines = 15; % Number of contour lines
vect_cc_xf = linspace(min(min(xfmat)),max(max(xfmat)),N_contour_lines);


 figure(1)
[xf_c,h_xfmat_c] = contourf(Vlin,nlin,xfmat,vect_cc_xf');
       clabel(xf_c,h_xfmat_c)
       colormap(flipud(colormap('gray')))
         grid on
         Title1 = strcat('Range (km) vs. V [m/s] & Engine revolutions [rps]');
         title(Title1)
         xlabel('V [m/s]')
         ylabel(' Engine revolutions [rps]')
         sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
         sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
         sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
         sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
        caxis([min(min(xfmat)) max(max(xfmat))]); %%%%%%%%%%%%%%%
         hold on 
         plot(Vlin,nconlin,'--r','LineWidth',2)
         legend('x_f (V,n)','T = D constrain')
         
figure(2)

mesh(Vlin,nlin,xfmat)


