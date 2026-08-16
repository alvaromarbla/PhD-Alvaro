function xi_dw = downwash_calc(alpha)

%% All should be in radians, not degrees

xi_0 = 0.01;
derxideralpha = 0.3;
xi_dw = xi_0 + derxideralpha * alpha;

end