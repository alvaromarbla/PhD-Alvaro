function [CL,CD] = aerodyn(alpha)
%This function provides the aerodynamic model of the ProVant-EMERGENTIa
%tilt-rotor aircraft.
p_CL = [0.043976 -0.004890 -0.2946314 0.377985 0.713766...
       -2.187085 -0.642635 3.345731 0.582];
p_CD = [0.0778734 -0.0586846 -0.485949 0.53177 1.049514...
       -2.376796 -0.852857 3.715949 0.0179749 0.08034];

CL = polyval(p_CL,alpha);
CD = polyval(p_CD,alpha);