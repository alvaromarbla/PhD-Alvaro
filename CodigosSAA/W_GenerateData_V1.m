function [W_data] = W_GenerateData_V1(D_data,Geo_data)
global SF

% Defines Estimation of Weights according to Cefiro III densities
%% identifies the aerodynamic surfaces being used
W1 = D_data.AC_CONFIGURATION.W1;
HTP = D_data.AC_CONFIGURATION.HTP;
VTP = D_data.AC_CONFIGURATION.VTP;
Can = D_data.AC_CONFIGURATION.Can;
Vee = D_data.AC_CONFIGURATION.Vee;

% Number of engines
n_eng = D_data.n_eng;

% Fuselage geometry
Surf_TOT = Geo_data.Surf_TOT;
f_f = Geo_data.f_f;
       
        
            % Weight of propulsion systems
            W_eng = 0.5; %
            m_engs = n_eng*W_eng;
            m_prop = 0.24;
            m_props = n_eng*m_prop;
            m_ESC = 0.145;
            m_ESCs = n_eng*m_ESC;
            m_prop = m_engs + m_props + m_ESCs;
            
            % subsystems
            n_servos = 6;
            m_servo = 0.024; % HITEC SERVOS
            m_servos = m_servo*n_servos;
            m_wiring = 0.200;
            m_metal = 0.250;
            m_subsystems = m_servos + m_wiring + m_metal;
            
            % ELECTRONIC ESTIMATION PROVANT 4.0 BRAZIL
            m_receiver = 0.016;
            m_NAVIO2 = 0.023;
            m_SONAR = 0.006;
            m_SHIELD1 = 0.150;
            m_SHIELD2 = 0.150;
            m_JETSON = 0.500;
            m_ACTUATOR_BUS = 0.050;
            m_POWER_SUPPLY = 0.100;
            m_GPS = 0.050;
            m_GPS_ANTENA = 0.173;
            m_IMU_ADIS = 0.048;
            m_PITOT= 0.025;
            m_wiring_electronics = 0.150;% Estimation
            m_systems = m_receiver + m_NAVIO2 + m_SONAR + m_SHIELD1 + m_SHIELD2 + m_JETSON + m_ACTUATOR_BUS + ...
                m_POWER_SUPPLY + m_GPS + m_GPS_ANTENA + m_IMU_ADIS + m_PITOT + m_wiring_electronics;
            
            m_landing_gear = 0;
            
            % Weight of PAX
            n_pax = 0;
            n_crew = 0;
            m_pax = 85; % mass for passanger
            m_crew = 85; % mass for crew
            m_luggage = 20; % mass for luggage and passanger
            m_pax_crew = (n_crew + n_pax)*m_luggage + n_crew*m_crew + n_pax*m_pax;
            
            % Weight of Payload
            m_cargo = 1*SF^3;
            m_payload = m_pax_crew + m_cargo;
            
            % Energy (fuel or batteries)
            m_batteries = 7;
            m_fuel = 0;
            m_energy = m_batteries + m_fuel;
            
            flag_landing_gear = 0; % determines if landing gear present (1 yes 0 no) 
        
        
   
        
        [rho_f,rho_fairing,rho_w,rho_HTP,rho_VTP,rho_tb]=CFIIIdensities; % calculated from Céfiro III
        rho_Vee = (rho_HTP + rho_VTP)/2;
        m_fus = rho_f*Surf_TOT ;            % Fuselage weight
        m_fairing = rho_fairing*Surf_TOT  ;             % Fairing weight
        rho_fus_fairing = (rho_f + rho_fairing)/2;
        rho_landing_gear = 0.04; % factor lineal landing gear
        rho_engine = 1; % factor lineal engine installation
        rho_misc = 0.10; % factor lineal msc
        
        if W1 == 1
            S_w1_s = Geo_data.S_w1_s;
            m_w1 = f_f*rho_w*S_w1_s  ;                           % Wing weight.
        else
            m_w1 = 0;
        end
        if HTP == 1
            S_w2_s = Geo_data.S_w2_s;
            m_HTP = f_f*rho_HTP*S_w2_s;                    %Peso HTP
        else
            m_HTP = 0;
        end
        if VTP == 1
            S_VTP_s = Geo_data.S_VTP_s;
            m_VTP = f_f*rho_VTP*S_VTP_s  ;                  %Peso HTP
        else
            m_VTP = 0;
        end
        if Can == 1
            S_w2_s = Geo_data.S_w2_s;
            m_Can = f_f*rho_w*S_w2_s;                  %Peso canard
        else
            m_Can = 0;
        end
        if Vee == 1
            S_w2_s = Geo_data.S_w2_s;
            m_Vee = f_f*rho_Vee*S_w2_s;                 %Peso cola en v
        else
            m_Vee = 0;
        end
        
        m_fus_fairing = f_f*rho_fus_fairing*Surf_TOT   ;                 %Peso fuselaje
        m_estructure = m_w1 + m_HTP + m_VTP + m_Can + m_Vee + m_fus_fairing;
        
        % Estimation
        m_prop = m_prop*rho_engine; % accounting for fairing, etc
        m_empty_estimation = m_estructure + m_prop + m_subsystems;
        MTOW_estimation = m_empty_estimation + m_payload + m_energy;
%         m_systems = f_f*rho_misc*MTOW_estimation;                       %Peso miscelaneos
        m_landing_gear = f_f*rho_landing_gear*MTOW_estimation*flag_landing_gear; 
        

        
        
        %Peso tren de aterrizaje
        m_empty = m_estructure + m_prop + m_landing_gear + m_subsystems + m_systems;
        m_TOW = m_empty + m_payload + m_energy;
        
        
        
%% Generate Output     
        
    
W_data.m_TOW = m_TOW;
W_data.m_energy = m_energy;
W_data.m_empty = m_empty;
W_data.m_payload = m_payload;
W_data.m_f_W0 = m_energy/m_TOW;
W_data.m_e_W0 = m_empty/m_TOW;
W_data.M_PAYLOAD.m_pax_crew = m_pax_crew;
W_data.M_PAYLOAD.m_cargo = m_cargo;
W_data.M_ENERGY.m_batteries = m_batteries;
W_data.M_ENERGY.m_fuel = m_fuel;
W_data.M_EMPTY.m_estructure = m_estructure;
W_data.M_EMPTY.m_prop = m_prop;
W_data.M_EMPTY.m_landing_gear = m_landing_gear;
W_data.M_EMPTY.m_subsystems = m_subsystems;
W_data.M_EMPTY.m_systems = m_systems;
W_data.M_ESTRUCTURE.m_w1 = m_w1;
W_data.M_ESTRUCTURE.m_HTP = m_HTP;
W_data.M_ESTRUCTURE.m_VTP = m_VTP;
W_data.M_ESTRUCTURE.m_Can = m_Can;
W_data.M_ESTRUCTURE.m_Vee = m_Vee;
W_data.M_ESTRUCTURE.m_fus_fairing = m_fus_fairing;
W_data.M_PROP.m_engs = m_engs;
W_data.M_PROP.m_ESCs = m_ESCs;
W_data.M_PROP.m_props = m_props;
W_data.M_SUBSYSTEMS.m_servos = m_servos;
W_data.M_SUBSYSTEMS.m_wiring = m_wiring;
W_data.M_SUBSYSTEMS.m_metal = m_metal;



end
