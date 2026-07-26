function slack = computeBoundSlack(X_hist, bounds)
    fields = fieldnames(bounds.vars);        % order matches how X was built
    n_iter = size(X_hist, 2);
    slack = struct();

    for k = 1:numel(fields)
        f  = fields{k};
        lb = bounds.vars.(f).min;
        ub = bounds.vars.(f).max;
        s_lo = X_hist(k,:) - lb;             % row k, same index used to build X
        s_hi = ub - X_hist(k,:);
        slack.(f).lower      = s_lo;
        slack.(f).upper      = s_hi;
        slack.(f).tightest   = min(s_lo, s_hi);
        slack.(f).normalized = min(s_lo, s_hi) ./ max(ub - lb, eps);
    end
end