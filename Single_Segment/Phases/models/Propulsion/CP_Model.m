function CP = CP_Model(V, n, alpha, epsilon, params)
    % Power coefficient model with zero-power constraint
    persistent a1_P b1_P c1_P d1_P e1_P b2_P c2_P d2_P b3_P c3_P b4_P;
    
    if isempty(a1_P)
        % Original coefficients from your code
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
    end
    
    % Calculate advance ratio and thrust angle
    J = V ./ (n .* params.prop.diameter);
    phi = alpha + epsilon;  % Consistent with original code
    
    % Polynomial expansion for CP
    CP = a1_P + b1_P.*J + c1_P.*J.^2 + d1_P.*J.^3 + e1_P.*J.^4 + ...
        abs(phi).*(b2_P.*J + c2_P.*J.^2 + d2_P.*J.^3) + ...
        abs(phi).^2..*(b3_P.*J + c3_P.*J.^2) + ...
        abs(phi).^3.*b4_P.*J;
    
end