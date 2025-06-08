% Power Model MDPI


a1_P = 0.0261518307541734;
b1_P = 0.0473735972985378;
c1_P = -0.16267474946046;
d1_P = 0.0247028469343899;
e1_P = 0.0306053713439883;

b2_P = -0.0762350484603968;
c2_P = 0.148580471912353;
d2_P = -0.0726017200715775;

b3_P = 0.0897273366920878;
c3_P = 0.0122602815262456;

b4_P = -0.0486029866039398;


CP = @(J,phi)  a1_P + b1_P*J + c1_P*J.^2 + d1_P*J.^3 + e1_P*J.^4 + ...
    phi.*(b2_P*J + c2_P*J.^2 + d2_P*J.^3) +  phi.^2.*(b3_P*J + c3_P*J.^2)...
    +  phi.^3.*(b4_P*J);

J_vec = linspace(0,1.2);
phi_vec = linspace(0,pi/2,180);


CP_mat = zeros(length(J_vec),length(phi_vec));
for ii = 1:length(J_vec)
    for jj = 1:length(phi_vec)
    CP_mat(ii,jj) = CP(J_vec(ii),phi_vec(jj));
    end
end


hFIG = figure;
fname = "CP_MDPI";

% Create contour
Levels = [-0.04 -0.03 -0.025 -0.02 -0.015 -0.01 -0.005 0 0.005 0.01 0.015 0.02 0.022 0.024 0.026  0.028  0.029 , 0.03, 0.031, 0.032, 0.033];
[c1,h1] = contourf(J_vec,phi_vec*180/pi,CP_mat',...
    'LevelList',Levels,'FaceAlpha',0.70);
colormap(flipud(colormap('gray')));
clabel(c1,h1,"Interpreter","latex",'FontSize',12);

hold on
% Create contour
contour(J_vec,phi_vec*180/pi,CP_mat','LineWidth',2,'LevelList',0,'LineColor',"r");

% Create ylabel
ylabel('Engine Tilt Angle $\phi$ [deg]');

% Create xlabel
xlabel('Advance Ratio J [-] ');

% Create title
title('Power Coefficient [-] vs. J [-] \& $\phi$ [deg]',"Interpreter","latex");

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