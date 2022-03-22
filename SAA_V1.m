%%% SENSITIVITY ANALYSIS ALGORITHM %%%


%% Requested changes

% CHANGE GLOBAL VARIABLES
% 

% LABELS: 

% Subsystem-Role-(...)-Version

% Subsystems: 

    % -> W : Weights
    % -> A : Aerodynamics
    % -> P : Propulsion
    % -> R : Performance
    % -> D : Design, Systems & Configuration
    % -> G : Graphics, Plots, and Result Processing 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global g     % Gravity
    global SF    % Scaling Factor
    
    g = 9.805;   %[m/s^2]
    SF = 1;      %[-]
    
    %% Select studies to perform
    horizontal_study = 0; % Set to 1 if you want to perform this study
    vertical_study   = 1; % Set to 1 if you want to perform this study
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate Aircraft Configuration & Geometry
    [D_data,Geo_data] = D_GenerateData_V1; 
    
    W_data = W_GenerateData_V1(D_data,Geo_data);
    
    P_data = P_GenerateData_V1(D_data);
    
    A_data = A_GenerateData_V1; %% Aerodynamics
    
    R_data = R_GenerateData_V1(W_data,Geo_data); %% Performance
    
    G_Options = G_Plot_Options_V1; 
        
   
   
    if horizontal_study == 1
    % Sensitivity Analysis Algorithm for Horizontal Flight
    [X_h,Y_h,Z_study_h,Fig_h] = P_SAA_Horizontal_V1(D_data, A_data,Geo_data,W_data,P_data,R_data,G_Options );
    end
    
    if vertical_study == 1
    % Sensitivity Analysis Algorithm for Vertical Flight
    [X_v,Y_v,Z_study_v,Fig_v] = P_SAA_Vertical_V1(D_data, A_data,Geo_data,W_data,P_data,R_data,G_Options );
    end
     