%function transition_corridor_thrust()
% This script calculates the transition corridor.

% Lower limit is calculated from condition of max operative angle of
% attack, solving n.



% Higher limit is calculated from max n and solving for alpha. 
params = Transition_Config();
S_ref = params.S_ref;
m     = params.mass;
rho   = params.rho;
g     = params.g;

eps_min = -5*pi/180;
eps_max = pi/2;
epsilon_vec = linspace(eps_min,eps_max,200);


%% Definition of aero forces
L = @(V,alpha) CL_Model(alpha) .* 0.5*rho*V.^2*S_ref;
D = @(V,alpha) CD_Model(alpha) .* 0.5*rho*V.^2*S_ref;

Eff_aero = @(V,alpha) L(V,alpha)./D(V,alpha);



%% Definition of functions to plot the corridor 
eps_fun = @(V,alpha) atan(Eff_aero(V,alpha).* (m*g./L(V,alpha) - 1)) - alpha;

T_fun   = @(V,alpha) (D(V,alpha).^2 +(m*g - L(V,alpha)).^2).^0.5;

Nlin = 200;
V_lin = linspace (0, 50, Nlin);
%% Calculate low limit 

alpha_lo = params.alpha_mp;


eps_lo = eps_fun(V_lin, alpha_lo);
T_lo   = T_fun (V_lin, alpha_lo);


%% Calculate high limit

alpha_hi = params.alpha_mi;

eps_hi = eps_fun(V_lin, alpha_hi);
T_hi   = T_fun (V_lin, alpha_hi);

figure(1)
plot(V_lin ,eps_lo *180/pi)
hold on
plot(V_lin, eps_hi*180/pi)

figure(2)
plot(V_lin ,T_lo )
hold on
plot(V_lin, T_hi)

%% Build closed corridor for fill
% Strip leading NaN (V=0 gives 0/0 in eps_fun)
mask_lo = ~isnan(eps_lo);
mask_hi = ~isnan(eps_hi);
V_lo_c = V_lin(mask_lo);   eps_lo_c = eps_lo(mask_lo);
V_hi_c = V_lin(mask_hi);   eps_hi_c = eps_hi(mask_hi);

% Keep only the portion of eps_lo that lies inside the corridor (eps >= -5 deg)
keep_lo  = eps_lo_c >= eps_min;
V_lo_c   = V_lo_c(keep_lo);
eps_lo_c = eps_lo_c(keep_lo);

% Keep only the portion of eps_hi that lies inside the corridor (eps >= -5 deg)
keep_hi  = eps_hi_c >= eps_min;
V_hi_c   = V_hi_c(keep_hi);
eps_hi_c = eps_hi_c(keep_hi);

% Closing segments
N_seg = 50;
% Bottom: short black segment at eps_min, between the two curves' floor endpoints
bot_seg_V = linspace(V_lo_c(end), V_hi_c(end), N_seg);
% Left: vertical at V=0, from eps_hi(V=0) up to eps_max (90 deg)
left_seg_eps = linspace(eps_hi_c(1), eps_max, N_seg);

%% Plot
hFIG = figure;

plot(V_lo_c, eps_lo_c*180/pi, 'LineWidth', 1.5, 'Color', 'r')
xlim([0 50])
ylim([-5 90])
hold on
plot(V_hi_c, eps_hi_c*180/pi, 'LineWidth', 1.5, 'Color', 'r')
plot(bot_seg_V, eps_min*ones(1,N_seg)*180/pi, 'LineWidth', 1.5, 'Color', 'k')
plot(zeros(1,N_seg), left_seg_eps*180/pi, 'LineWidth', 1.5, 'Color', 'k')

% Polygon walk (counter-clockwise):
%   eps_lo (left -> right, ending on the floor)
%   bottom segment along eps_min (V_lo_c(end) -> V_hi_c(end))
%   eps_hi reversed (right -> left, ending at V=0)
%   left segment at V=0 (eps_hi_c(1) -> eps_max), then implicitly closes to eps_lo_c(1)
X_fill = [V_lo_c, bot_seg_V, fliplr(V_hi_c), zeros(1,N_seg)];
Y_fill = [eps_lo_c*180/pi, eps_min*ones(1,N_seg)*180/pi, ...
          fliplr(eps_hi_c)*180/pi, left_seg_eps*180/pi];

ylabel('$\varepsilon$ [deg]');
xlabel('V [m/s]');
grid on
title('Transition corridor (thrust formulation)');

fill(X_fill, Y_fill, [0.7,0.7,0.7], 'FaceAlpha', .2)

fname = "transition_corridor_thrust";

picturewidth = 20;
hw_ratio = 0.65;
set(findall(hFIG,'-property','FontSize'),'FontSize',16)
text(5,30,"Low Speed Limit",'FontSize',12);
text(35,75,"High Speed Limit",'FontSize',12);
set(findall(hFIG,'-property','Box'),'Box','off')
set(findall(hFIG,'-property','Interpreter'),'Interpreter','latex')
set(findall(hFIG,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hFIG,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hFIG,'Position');
set(hFIG,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(hFIG,fname,'-dpdf','-vector','-fillpage')