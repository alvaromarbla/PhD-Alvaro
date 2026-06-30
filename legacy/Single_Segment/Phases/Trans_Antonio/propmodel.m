function [C_T,C_P] = propmodel(J,phi)
%This function provides the thrust coefficient and the power coefficient of
%a fixed pitch propeller for given values of the advance ratio (J) and the
%aerodynamic tilt angle (phi). J and phi must be equal-size vectors.

C_T = 0.0735 - 0.0311*J - 0.2497*J.^2 + 0.1430*J.^3 + 0.0261*J.^4 + ...
    -0.0982*phi.*J + 0.20127*phi.*J.^2 - 0.1737*phi.*J.^3 + ...
    0.1562*phi.^2.*J + 0.0368*phi.^2.*J.^2 - 0.0478*phi.^3.*J;

C_P = 0.0261 + 0.0473*J - 0.1626*J.^2 + 0.0247*J.^3 + 0.0306*J.^4 + ...
    -0.0762*phi.*J + 0.1485*phi.*J.^2 - 0.0726*phi.*J.^3 + ...
    0.0897*phi.^2.*J + 0.0122*phi.^2.*J.^2 - 0.0486*phi.^3.*J;

%%%%%%%%%%%
% J_vec = (0:1e-3:1)';
% phi_vec = pi/2*(0:1e-3:1)';
% 
% CT_mat = zeros(length(J_vec),length(phi_vec));
% CP_mat = CT_mat;
% 
% for i=1:length(J_vec)
%     for j=1:length(phi_vec)
%         [C_Tij,C_Pij] = propmodel(J_vec(i),phi_vec(j));
%         CT_mat(i,j) = C_Tij;
%         CP_mat(i,j) = C_Tij;
%     end
% end
% 
% figure, [C,h] = contour(J_vec,phi_vec,CT_mat',(-0.04:0.02:0.14)'); clabel(C,h)