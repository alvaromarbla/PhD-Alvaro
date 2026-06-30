function [c, ceq] = Power_Constraint(CP)

%% NOT USED 
    % Enforces C_P = 0 (no power consumption)
    ceq = CP;  % Equality constraint CP = 0
    c = [];     % No inequality constraints
    
end