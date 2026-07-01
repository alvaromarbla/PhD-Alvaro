function slack = computeBoundSlack(X_hist, bounds)
    idx = struct('V',1,'gamma',2,'alpha',3,'epsilon',4,'n',5);
    fields = fieldnames(idx);
    n_iter = size(X_hist, 2);
    slack = struct();
    for k = 1:numel(fields)
        f = fields{k};
        i = idx.(f);
        lb = bounds.(f).min;
        ub = bounds.(f).max;
        s_lo = X_hist(i,:) - lb;
        s_hi = ub - X_hist(i,:);
        slack.(f).lower = s_lo;
        slack.(f).upper = s_hi;
        slack.(f).tightest = min(s_lo, s_hi);      % distance to nearest bound
        slack.(f).normalized = min(s_lo, s_hi) ./ max(ub - lb, eps); % 0 = at bound, 0.5 = centered
    end
end