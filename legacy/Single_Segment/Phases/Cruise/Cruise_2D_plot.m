function fig = Cruise_2D_plot(sol,fval,params)
%% Generate numerical Matrix
Nmat = 10;
Vlin = linspace(17,35,Nmat);
nlin = linspace(35,95,Nmat);

% Preallocate matrix for speed
xfmat = zeros(Nmat,Nmat);
% Define x as input for Cruise_objective
x(2) = sol(2); %alpha. Not needed for cruise
x(3) = 0; %epsilon
% Loop in speed
for  ii = 1: Nmat
    % Loop in revolutions
    for jj = 1:Nmat
        x(1) = Vlin(jj);
        x(4) = nlin(ii);
        xfmat(ii,jj) = -Cruise_Objective(x, params)*1e-03;
    end
end
xfmat
%% Extract solution

V_sol = sol(1);
n_sol = sol(4);
x_sol = -fval*1e-03; % To make into km

%% Plot contours
N_contour_lines = 12; % Number of contour lines
%vect_cc_xf = [5,10,15,20,25,30,40,50,60,70,80];

% The usual command for vect_cc_xf is below, but we use the line above
% because we know the result (for prettier result!)
vect_cc_xf = linspace( 0,200 ,N_contour_lines);
%min(min(xfmat)),max(max(xfmat))

fig = figure(1);
[xf_c,h_xfmat_c] = contourf(nlin,Vlin,xfmat,vect_cc_xf');
       clabel(xf_c,h_xfmat_c)
       colormap(flipud(colormap('gray')))
         grid on
         Title1 = strcat(['Range [km] vs. V [m/s] & Engine rev. [rps]. Max Range is ', num2str(x_sol), ' km']); 
         title(Title1)
         xlabel('Engine revolutions [rps] ')
         ylabel('V [m/s]')
         sidebar = colorbar;  %%%%%%%%%%%%%%%%%%%
         sidebar.Label.FontSize = 11;   %%%%%%%%%%%%%%%%%%%
         sidebar.Label.FontWeight = 'bold';   %%%%%%%%%%%%%%%%%%%
         sidebar.Label.Position(1) = 2;   %%%%%%%%%%%%%%%%%%%
         caxis([min(min(xfmat)) max(max(xfmat))]); %%%%%%%%%%%%%%%
         hold on 
         %plot(nconlin,Vlin,'--r','LineWidth',2)
         %plot(nconlin_full,Vlin,'-.m','LineWidth',2)
         %plot(nTmaxlin,Vlin,'-b','LineWidth',2.5)
         %yline(Vstall,'-.b','LineWidth',2)
         xline(params.prop.max_rps,':b','LineWidth',2)
         %plot(nsol_A,Vsol,'*g','LineWidth',2)
         plot(n_sol,V_sol,'og','LineWidth',2.5)
         legend('x_f (V,n)','T = D constrain. Full Model','Min. Operative Speed','Max RPS','Optimal solution. Full Model','Location','northwest')
