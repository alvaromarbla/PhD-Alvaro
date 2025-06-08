function P =  prop_CP_1365(V, n , phi, params)

    persistent a1_P b1_P c1_P d1_P e1_P b2_P c2_P d2_P b3_P c3_P b4_P J_min J_max phi_min phi_max ;
    
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
    
        J_min = 0;
        J_max = 1.2; 

        phi_min = 0;
        phi_max = pi/2;

    end
    
    %% Warnings on validity of the model

    % Calculate advance ratio and thrust angle
    J = V ./ (n .* params.prop.diameter);
    
    if any(J)<J_min || any(J)>J_max
         warning('Advance ratio J = %.3f is outside the model''s valid range [%.2f, %.2f]. Results may not be valid.', J, J_min, J_max);
    end 

    if any(abs(phi))<phi_min || any(abs(phi))>phi_max
        warning('Incidence angle phi = %.3f rad is outside the model''s valid range [%.2f, %.2f] rad. Results may not be valid.', abs(phi), phi_min, phi_max);
    end 

    % Polynomial expansion for CP
    CP = a1_P + b1_P.*J + c1_P.*J.^2 + d1_P.*J.^3 + e1_P.*J.^4 + ...
        abs(phi).*(b2_P.*J + c2_P.*J.^2 + d2_P.*J.^3) + ...
        abs(phi).^2..*(b3_P.*J + c3_P.*J.^2) + ...
        abs(phi).^3.*b4_P.*J;

    % Abs value in phi ensures symmetry w.r.t. phi if it's negative

    % Convert to Power [W]
    P = CP .* params.prop.num_engines .* params.rho .* n.^3 .* params.prop.diameter^5;

end