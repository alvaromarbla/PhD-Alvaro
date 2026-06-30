function CL = CL_Model(alpha)
    % 8th-order polynomial for C_L(α) from wind tunnel data
    p_CL_ac = [0.582; 3.345731; -0.642635; -2.187085; 0.713766; ...
               0.377985; -0.2946314; -0.004890; 0.043976];
    
    % Manual polynomial evaluation (matches original code structure)
    CL = p_CL_ac(9)*alpha.^8 + p_CL_ac(8)*alpha.^7 + p_CL_ac(7)*alpha.^6 + ...
         p_CL_ac(6)*alpha.^5 + p_CL_ac(5)*alpha.^4 + p_CL_ac(4)*alpha.^3 + ...
         p_CL_ac(3)*alpha.^2 + p_CL_ac(2)*alpha + p_CL_ac(1);
end