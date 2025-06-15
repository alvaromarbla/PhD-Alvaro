function [initial_cond] = cruise_initial_cond (previous_cond, params)


x0.V = 24; %  [m/s]
x0.gamma = deg2rad(0); %[rad]
x0.alpha = deg2rad(5); %[rad]
x0.epsilon = deg2rad(5); %[rad]
x0.n = 30; %[rps]

% X(1) = V
% X(2) = gamma
% X(3) = alpha
% X(4) = epsilon
% X(5) = n

initial_cond = [x0.V, x0.gamma, x0.alpha, x0.epsilon, x0.n];

end 