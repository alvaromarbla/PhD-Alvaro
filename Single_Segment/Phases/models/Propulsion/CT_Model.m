function [T, phi] = CT_Model(V, n, alpha, epsilon, params)
    % Thrust coefficient model with thrust vector angle
    persistent a1_T b1_T c1_T d1_T e1_T b2_T c2_T d2_T b3_T c3_T b4_T;
    
    if isempty(a1_T)
        % Original coefficients from your code
        a1_T = 0.0735880531010883;
        b1_T = -0.0311758018412727;
        c1_T = -0.249744726429543;
        d1_T = 0.143084420694372;
        e1_T = 0.0261032283758581;
        b2_T = -0.0982459868751664;
        c2_T = 0.20127470719351;
        d2_T = -0.173738749783189;
        b3_T = 0.156239779501715;
        c3_T = 0.0368592048084175;
        b4_T = -0.0478709034281346;
    end
    
    % Calculate advance ratio and thrust angle
    J = V./ (n .* params.prop.diameter);
    phi = alpha + epsilon;  % Original code's definition
    
    % Polynomial expansion for CT
    CT = a1_T + b1_T.*J + c1_T.*J.^2 + d1_T.*J.^3 + e1_T.*J.^4 + ...
        abs(phi).*(b2_T.*J + c2_T.*J.^2 + d2_T.*J.^3) + ...
        abs(phi).^2..*(b3_T.*J + c3_T.*J.^2) + ...
        abs(phi).^3..*b4_T.*J;
    
    % Convert to physical thrust [N]
    T = CT .* params.prop.num_engines .* params.rho .* n.^2 .* params.prop.diameter^4;
end