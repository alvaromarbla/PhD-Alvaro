% Thrust Model MDPI Paper

a1_T = 0.0735880531010883;
b1_T = -0.0311758018412727;
c1_T = -0.249744726429543;
d1_T = 0.143084420694372;
e1_T = 0.0261032283758581;

b2_T = -0.0982459868751664;
c2_T = 0.20127470719351;
d2_T = -0.173738749783189;

b3_T = 0.156239779501715;
c3_T = 0.0368592048084175;

b4_T = -0.0478709034281346;


CT = @(J,phi)  a1_T + b1_T*J + c1_T*J.^2 + d1_T*J.^3 + e1_T*J.^4 + ...
    phi.*(b2_T*J + c2_T*J.^2 + d2_T*J.^3) +  phi.^2.*(b3_T*J + c3_T*J.^2)...
    +  phi.^3.*(b4_T*J);

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

J_vec = linspace(0,2);
phi_vec = linspace(0,pi/2,180);


eta_mat = zeros(length(J_vec),length(phi_vec));
for ii = 1:length(J_vec)
    for jj = 1:length(phi_vec)
    eta_mat(ii,jj) = cos(phi_vec(jj)).*J_vec(ii).*CT(J_vec(ii),phi_vec(jj))./CP(J_vec(ii),phi_vec(jj));
    end
end

hFIG = figure;
fname = "eta_MDPI";

% Create contour
Levels = [-1 -0.5 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.3 1.4 1.5 2 3 4 5 6 ];
[c1,h1] = contourf(J_vec,phi_vec*180/pi,eta_mat',...
    'LevelList',Levels,'FaceAlpha',0.70);
colormap(flipud(colormap('gray')));
clabel(c1,h1,"Interpreter","latex",'FontSize',12);

hold on
% Create contour
% contour(J_vec,phi_vec*180/pi,eta_mat','LineWidth',2,'LevelList',0,'LineColor',"r");

% Create ylabel
ylabel('Engine Tilt Angle $\phi$ [deg]');

% Create xlabel
xlabel('Advance Ratio J [-] ');

% Create title
title('Efficiency [-] vs. J [-] \& $\phi$ [deg]',"Interpreter","latex");

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