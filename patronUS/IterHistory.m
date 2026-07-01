classdef IterHistory < handle
    properties
        X_hist = []
        g_hist = []
        count = 0
    end
    methods
        function record(obj, Xval, gval)
            obj.count = obj.count + 1;
            obj.X_hist(:, obj.count) = Xval;
            obj.g_hist(:, obj.count) = gval;
        end
    end
end