function CD = CD_wholeAC_1wing(alpha)
    % 9th-order polynomial for C_D(α) from wind tunnel data
    p_CD_ac = [0.08034; 0.0179749; 3.715949; -0.852857; -2.376796; ...
               1.049514; 0.53177; -0.485949; -0.0586846; 0.0778734];
    
    % Manual polynomial evaluation preserving original order
    CD = p_CD_ac(10)*alpha.^9 + p_CD_ac(9)*alpha.^8 + p_CD_ac(8)*alpha.^7 + ...
         p_CD_ac(7)*alpha.^6 + p_CD_ac(6)*alpha.^5 + p_CD_ac(5)*alpha.^4 + ...
         p_CD_ac(4)*alpha.^3 + p_CD_ac(3)*alpha.^2 + p_CD_ac(2)*alpha + p_CD_ac(1);
end