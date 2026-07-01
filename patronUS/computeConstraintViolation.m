function viol = computeConstraintViolation(g_hist, n_ineq, n_eq)
    c_hist   = g_hist(1:n_ineq, :);
    ceq_hist = g_hist(n_ineq+1:n_ineq+n_eq, :);

    viol.ineq_per_row = max(0, c_hist);        % violation per inequality, per iter
    viol.eq_per_row   = abs(ceq_hist);         % violation per equality, per iter

    viol.ineq_max = max(viol.ineq_per_row, [], 1);   % worst inequality violation per iter
    viol.eq_max   = max(viol.eq_per_row, [], 1);     % worst equality violation per iter
    viol.total    = max(viol.ineq_max, viol.eq_max); % analogous to IPOPT's inf_pr
end