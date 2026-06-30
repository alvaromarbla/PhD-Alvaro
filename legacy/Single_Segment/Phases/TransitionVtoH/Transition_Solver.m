function [results, total_energy] = Transition_Solver(params)
    % Phase 1: Fixed α, varying ε
    %phase1_params = params;
    % phase1_params.phi_0 = pi/2;  % Initial vertical tilt
    % t_est = params.phi_0/params.epsilon_slope;
    
    alpha_mi = params.alpha_mi;
    alpha_t  = params.alpha_t;
    V_f      = params.V_f;
    tau_tr   = params.tau_tr;
    options = odeset('RelTol', params.RelTol,...
                    'Events', @(t,y) Phase1_Events(t,y,params));
    [t_trans, y_trans] = DAESolver(@(t,y) Phase1_DAE(t,y,params),...
                        params.tspan, [0; 0; 0; params.prop.n_hov], options);
    
    % Combine results
    results.t = t_trans;
    results.x = y_trans(:,1);
    results.V = y_trans(:,2);
    results.E = y_trans(:,3);
    results.n = y_trans(:,4);
    results.epsilon_vec = pi/2*(1-t_trans/tau_tr).^2;
    results.alpha_vec = alpha_t + (alpha_mi-alpha_t)*exp(-5*t_trans/tau_tr);
    results.params = params;
    
    results.L =  CL_Model(results.alpha_vec) .* 0.5*params.rho.*results.V.^2.*params.S_ref;
    results.D =  CD_Model(results.alpha_vec) .* 0.5*params.rho.*results.V.^2.*params.S_ref;

    [results.T, ~] = CT_Model(results.V, results.n, results.alpha_vec, results.epsilon_vec, params);
    results.P = CP_Model(results.V, results.n, results.alpha_vec, results.epsilon_vec, params)*params.rho.*params.prop.num_engines.*results.n.^3.*params.prop.diameter^5;

    % results.epsilon_vec = 0.1*pi/2*exp(-4*t_trans) + ...
    % 0.9*pi/2*(1 - exp(min(results.V, V_f)/4) ./ exp(V_f/4));
    % results.alpha_vec   = alpha_mi*(1-exp(-4*t_trans)) + (alpha_t-alpha_mi)*...
    % (exp(min(results.V, V_f)/4 )./ exp(V_f/4));


    % Calculate total energy consumption
    total_energy = results.E(end) - results.E(1);

    % Plot transition corridor and solution

    if params.flag_corridor ==1
        %transition_corridor()
        %hold on

        hFIG1 = figure(7);
        fname = "distance_trans";
        plot(results.t,results.x,"Color", "b", "LineWidth",1.5)
        xlabel("t [s]")
        ylabel("x [m]")
        title("Distance [m] vs Time [s] evolution. ")
        grid on
        picturewidth = 20;
        hw_ratio = 0.65;
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex')
        set(findall(hFIG1,'-property','FontSize'),'FontSize',16)
        set(findall(hFIG1,'-property','Box'),'Box','off') % optional
        set(findall(hFIG1,'-property','Interpreter'),'Interpreter','latex') 
        set(findall(hFIG1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
        set(hFIG1,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
        pos = get(hFIG1,'Position');
        set(hFIG1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        print(hFIG1,fname,'-dpdf','-vector','-fillpage')


        hFIG2 = figure(2);
        fname = "speed_trans";
        plot(results.t,results.V,"Color", "b", "LineWidth",1.5)
        xlabel("t [s]")
        ylabel("V [m/s]")
        hold on
        grid on
        yline(V_f,"--r", "LineWidth",1.5 )
        title("Speed [m/s] vs Time [s] evolution. ")
        legend("Evolution of the speed","$V_f$",'Location','east')
        picturewidth = 20;
        hw_ratio = 0.65;
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex')
        set(findall(hFIG2,'-property','FontSize'),'FontSize',16)
        set(findall(hFIG2,'-property','Box'),'Box','off') % optional
        set(findall(hFIG2,'-property','Interpreter'),'Interpreter','latex') 
        set(findall(hFIG2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
        set(hFIG2,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
        pos = get(hFIG2,'Position');
        set(hFIG2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        print(hFIG2,fname,'-dpdf','-vector','-fillpage')


        hFIG3 = figure(3);
        fname = "n_trans";
        plot(results.t,results.n,"Color", "b", "LineWidth",1.5)
        xlabel("t [s]")
        ylabel("n [rev/s]")
        hold on
        grid on
        yline(params.prop.rps_max,"--r","LineWidth",1.5 )
        title("Engine running speed [rad/s] vs Time [s] evolution. ")
        legend("Evolution of the engine running speed","$n_{max}$", "Location" , "west")
        picturewidth = 20;
        hw_ratio = 0.65;
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex')
        set(findall(hFIG3,'-property','FontSize'),'FontSize',16)
        set(findall(hFIG3,'-property','Box'),'Box','off') % optional
        set(findall(hFIG3,'-property','Interpreter'),'Interpreter','latex') 
        set(findall(hFIG3,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
        set(hFIG3,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
        pos = get(hFIG3,'Position');
        set(hFIG3,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        print(hFIG3,fname,'-dpdf','-vector','-fillpage')    


        hFIG4 = figure(4);
        fname = "force_trans"; 
        plot(results.t,results.T, "b", "LineWidth",1.5 )
        xlabel("t [s]")
        ylabel("Forces [N]")
        hold on 
        grid on
        plot(results.t,results.L, "--r", "LineWidth",1.5 )
        plot(results.t,results.D, "-.k", "LineWidth",1.5)
        title("Main acting forces [N] vs Time [s] evolution")
        legend("Engine Thrust", "Aerodynamic Lift", "Aerodynamic Drag", "Location" , "East")
        picturewidth = 20;
        hw_ratio = 0.65;
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex')
        set(findall(hFIG4,'-property','FontSize'),'FontSize',16)
        set(findall(hFIG4,'-property','Box'),'Box','off') % optional
        set(findall(hFIG4,'-property','Interpreter'),'Interpreter','latex') 
        set(findall(hFIG4,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
        set(hFIG4,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
        pos = get(hFIG4,'Position');
        set(hFIG4,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        print(hFIG4,fname,'-dpdf','-vector','-fillpage') 



        hFIG5 = figure(5);
        fname = "consumption_trans";
        plot(results.t,results.E/params.E_batt*100,"b", "LineWidth",1.5 )
        xlabel("t [s]")
        ylabel("$E_{consumed}$ [\%]")
        hold on
        grid on
        yyaxis right
        axis([0, 12, 0, 3])
        ylabel("$P_{consumed}$ [kW]")
        plot(t_trans,results.P*1e-03, "--k", "LineWidth",1.5 )
        title("Energy [\%] and Power [kW] consumption vs Time [s] evolution" )
        legend("Percentage of energy consumption w.r.t. total", "Power consumption","Location", "northwest")
        picturewidth = 20;
        hw_ratio = 0.65;
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex')
        set(findall(hFIG5,'-property','FontSize'),'FontSize',16)
        set(findall(hFIG5,'-property','Box'),'Box','off') % optional
        set(findall(hFIG5,'-property','Interpreter'),'Interpreter','latex') 
        set(findall(hFIG5,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
        set(hFIG5,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
        pos = get(hFIG5,'Position');
        set(hFIG5,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        print(hFIG5,fname,'-dpdf','-vector','-fillpage') 


        hFIG6 = figure(6);
        fname = "angles_trans";
        plot(results.t,results.epsilon_vec*180/pi,"b", "LineWidth",1.5 )
        xlabel("t [s]")
        ylabel("Angles [deg]")
        hold on
        grid on
        plot(t_trans,results.alpha_vec*180/pi,"--r", "LineWidth",1.5 )
        plot(t_trans,(results.alpha_vec+results.epsilon_vec)*180/pi,"-.k", "LineWidth",1.5 )
        title("Main Angles [deg] in action vs Time [s]")
        legend("Engine tilt angle $\varepsilon$", "Angle of attack $\alpha$", "Engine inclination w.r.t. airflow $\phi$",'Location','northeast')
        picturewidth = 20;
        hw_ratio = 0.65;
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex')
        set(findall(hFIG6,'-property','FontSize'),'FontSize',16)
        set(findall(hFIG6,'-property','Box'),'Box','off') % optional
        set(findall(hFIG6,'-property','Interpreter'),'Interpreter','latex') 
        set(findall(hFIG6,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
        set(hFIG6,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
        pos = get(hFIG6,'Position');
        set(hFIG6,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        print(hFIG6,fname,'-dpdf','-vector','-fillpage') 

       

    

    end



end


% Phase1_Events.m
function [value, isterminal, direction] = Phase1_Events(t, y, params)
    value = y(2) - params.V_f; %params.V_min_ope;  % Reaches operational speed
    isterminal = 1;
    direction = 1;
end