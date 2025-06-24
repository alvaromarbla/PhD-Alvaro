function matrix_conditions = cruise_initial_cond (lb, ub, N_samples)
% CRUISE_INITIAL_COND   Generate Latin-Hypercube initial conditions
%
%
%   Outputs an [nVars × N_samples] matrix where:
%     • each row i is samples for variable i,
%     • each column j is the full-variable vector for sample j.
%
%   Inputs:
%     lb, ub                 - vectors of lower and upper bounds
%     N_samples              - number of LHS samples

% Identify which dims actually vary

varying    = (lb ~= ub);
nVary      = sum(varying);

% Generate LHS in [0,1] for only the varying dims
L = lhsdesign(N_samples, nVary);

% Preallocate: rows = variables, cols = samples
nVars = numel(lb);
matrix_conditions = zeros(nVars, N_samples);


% Vectorized scaling:
%   For each varying i: row varIdx(i) = L(:,i)' * (ub-lb) + lb
diffs = (ub(varying) - lb(varying));
matrix_conditions(varying, :) = bsxfun(@times, L', diffs) + lb(varying);

% Fill fixed-variable rows with their constant bound
fixedIdx = find(~varying);
for idx = fixedIdx.'
    matrix_conditions(idx, :) = lb(idx);
end

% x0.V = 24; %  [m/s]
% x0.gamma = deg2rad(5); %[rad]
% x0.alpha = deg2rad(5); %[rad]
% x0.epsilon = deg2rad(50); %[rad]
% x0.n = 30; %[rps]

% X(1) = V
% X(2) = gamma
% X(3) = alpha
% X(4) = epsilon
% X(5) = n

% matrix_conditions = [x0.V, x0.gamma, x0.alpha, x0.epsilon, x0.n]';

end