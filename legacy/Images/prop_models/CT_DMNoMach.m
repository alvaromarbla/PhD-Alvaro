% Thrust Model Diego Manzano No Mach

a1_T = 0.0746411310089907 ;
b1_T = -0.0511217219165455;
c1_T = -0.185616842464272;
d1_T = 0.099588704442492;
e1_T = -0.00475509466498915;

a2_T = 0.00628324456145052;
b2_T = -0.0641445035520659;
c2_T = 0.0898897428397212;
d2_T = -0.0684236614384807;

a3_T = -0.0346992361203213;
b3_T = 0.165130565299384;
c3_T = 0.0339434014319326;

a4_T = 0.0402868761263794;
b4_T = -0.0545228637854472;

a5_T = -0.0142867202072431;


CT = @(J,phi)  a1_T + b1_T*J + c1_T*J.^2 + d1_T*J.^3 + e1_T*J.^4 + ...
    phi.*(a2_T + b2_T*J + c2_T*J.^2 + d2_T*J.^3) +  phi.^2.*(a3_T + b3_T*J + c3_T*J.^2)...
    +  phi.^3.*(a4_T + b4_T*J) + phi.^4.*a5_T;

J_vec = linspace(0,1);
phi_vec = linspace(0,pi/2,180);


CT_mat = zeros(length(J_vec),length(phi_vec));
for ii = 1:length(J_vec)
    for jj = 1:length(phi_vec)
    CT_mat(ii,jj) = CT(J_vec(ii),phi_vec(jj));
    end
end


hFIG = figure;
fname = "CT_DMNoMach";

% Create contour
[c1,h1] = contourf(J_vec,phi_vec*180/pi,CT_mat',...
    'LevelList',[-0.04 -0.02 0 0.02 0.04 0.06 0.07 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28]);
clabel(c1,h1,"Interpreter","latex");

hold on
% Create contour
contour(J_vec,phi_vec*180/pi,CT_mat','LineWidth',2,'LevelList',0,'LineColor',"r");

% Create ylabel
ylabel('Engine Tilt Angle $\phi$ [deg]');

% Create xlabel
xlabel('Advance Ratio J [-] ');

% Create title
title('CT Model Diego M. No Mach [-] vs. J [-] \& $\phi$ [deg]',"Interpreter","latex");

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