% Power Model 1410


a1_P = 0.0334878516820497;
b1_P = 0.0462334015755113;
c1_P = -0.0490708225621856;
d1_P = -0.135072940509969;
e1_P = 0.0703737279878982;

a2_P = 0.00477674951107802;
b2_P = -0.119198383260838;
c2_P = 0.159426123331315;
d2_P = -0.0147991437397226;

a3_P = 0.0216950039616885;
b3_P = 0.115262163726395;
c3_P = -0.0308423257753396;

a4_P = -0.0406171047062065;
b4_P = -0.0428659151604016;

a5_P = 0.0168435067341552;


CP = @(J,phi)  a1_P + b1_P*J + c1_P*J.^2 + d1_P*J.^3 + e1_P*J.^4 + ...
    phi.*(a2_P + b2_P*J + c2_P*J.^2 + d2_P*J.^3) +  phi.^2.*(a3_P + b3_P*J + c3_P*J.^2)...
    +  phi.^3.*(a4_P + b4_P*J) + phi.^4.*a5_P;

J_vec = linspace(0,1);
phi_vec = linspace(0,pi/2,180);


CP_mat = zeros(length(J_vec),length(phi_vec));
for ii = 1:length(J_vec)
    for jj = 1:length(phi_vec)
    CP_mat(ii,jj) = CP(J_vec(ii),phi_vec(jj));
    end
end


hFIG = figure;
fname = "CP_1410";

% Create contour
[c1,h1] = contourf(J_vec,phi_vec*180/pi,CP_mat',...
    'LevelList',[-0.04 -0.03 -0.025 -0.02 -0.015 -0.01 -0.005 0 0.005 0.01 0.015 0.02 0.025 0.03 0.032 0.034 0.036 0.038 0.04 0.045 0.050 ]);
clabel(c1,h1,"Interpreter","latex");

hold on
% Create contour
contour(J_vec,phi_vec*180/pi,CP_mat','LineWidth',2,'LevelList',0,'LineColor',"r");

% Create ylabel
ylabel('Engine Tilt Angle $\phi$ [deg]');

% Create xlabel
xlabel('Advance Ratio J [-] ');

% Create title
title('CP Model 1410 [-] vs. J [-] \& $\phi$ [deg]',"Interpreter","latex");

% box(axes,'on');
% grid(axes,'on');
% axis(axes,'tight');
% hold(axes,'off');
% % Set the remaining axes properties
% set(axes,'BoxStyle','full','CLim',[-0.04 0.34],'Layer','top');
% % Create colorbar
colorbar;
% 
picturewidth = 20;
hw_ratio = 0.65;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex')
set(findall(hFIG,'-property','FontSize'),'FontSize',16)
set(findall(hFIG,'-property','Box'),'Box','off') % optional
set(findall(hFIG,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hFIG,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hFIG,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hFIG,'Position');
set(hFIG,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hFIG,fname,'-dpdf','-vector','-fillpage')