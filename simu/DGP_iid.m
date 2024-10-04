function [ Y, Y_test, Theta_tt, breaks ] = DGP_iid(T, p, m, prob, tau, Nsim, regimes)
%% Generate iid data set.
%
% [Y, Y_test, Theta, breaks] = DGP_iid(T, p, m, prob, regimes);
%
% - Input:
%   @T:         Sample size.
%   @p:         Dimension.
%   @m:         Number of break points.
%   @prob:      The sparsity of each true Theata.
%   @regimes:   'costant' or 'non-constant' to control whether the length of
%               each regime is same.
%
% - Output:
%   @Y:         Training set of size T x p.
%   @Y_test:    Test set of size T x p x 4.
%   @Theta:     True Theta.
%   @breaks:    True breaks.


arguments (Input)
    T (1, 1) double = 100 % sample size
    p (1, 1) double = 10 % dimension
    m (1, 1) double = 2 % number of break points
    prob double = 0.8 % the sparsity of each true Theta's.
    tau (1, 1) double = 14 % the rescaling factor of Theta.
    Nsim (1, 1) int32 = 10 % the number of simulated validation datasets.
    regimes string = "non-constant" % "non-constant" or "constant".
end

assert(any(strcmp(regimes, {'constant', 'non-constant'})), ...
    "The 'regimes' must be 'constant' or 'non-constant'.");

% Nblocks: number of distinct blocks of data
Nblocks = m+1;
eps = 1 / (m + 8);

% The length of the input prob.
len_prob = length(prob);

% If the prob is a scalar, i.e., the sparsities for all regimes are the
% same, then copy it for Nblocks times.
if len_prob == 1
    prob = prob * ones(Nblocks, 1); % sparsity = round((1-prob)*p*(p-1)/2);
else
    % Otherwise, if the length is not one but less then Nblocks, then fill
    % the remaining part by 0.8.
    if len_prob < Nblocks
        prob = [reshape(prob, [], 1); 0.8 * ones(Nblocks - len_prob)]; % sparsity = round((1-prob)*p(p-1)/2);
    end
end

switch regimes
    case 'non-constant'
        % specify the break dates such that the length of each block is
        % at least eps*n
        cond = true; i = 0;
        while cond
            i = i+1;
            breaks = sort(randi([round(eps*T)+1 T-round(eps*T)],m,1));
            diff = breaks(2:end)-breaks(1:end-1);
            cond = (min(diff)<eps*T);
            if i > 10^5
                break
            end
        end
    case 'constant'
        % specify the break dates such that the sample size for each regime
        % is similar
        len = round(T/(m+1)); breaks = ((1:m)*len)';
end

Theta = zeros(p,p,m+1);
a1 = -1; a2 = 1; b1 = 0.1; b2 = 0.5;
for k = 1:m+1
    sparsity = round(p*(p-1)/2*prob(k));
    Theta(:,:,k) = simulate_sparse_precision(p,sparsity,prob(k),a1,a2,b1,b2); 
end
Theta = Theta/tau;
% Generate the data based on a Gaussian distribution with
% variance-covariance matrix Sigma = inv(Theta)
% Generate the data based on a Gaussian distribution with
% variance-covariance matrix Sigma = inv(Theta)
if m>0 % if one break at least
    Y = zeros(T, p);
    ind_ = [0; breaks; T];
    for i = 1:m + 1
        Y(ind_(i) + 1 : ind_(i + 1), :) = mvnrnd(zeros(p, 1), ...
            inv(Theta(:, :, i)), ind_(i + 1) - ind_(i));
    end
    
    Y_test = zeros(T, p, Nsim);
    for tt = 1 : Nsim
        for i = 1 : m + 1
            Y_test(ind_(i) + 1 : ind_(i + 1), :, tt) = mvnrnd(zeros(p, 1), ...
                inv(Theta(:, :, i)), ind_(i + 1) - ind_(i));
        end
    end
    
    Theta_tt = zeros(p, p, T);
    Theta_tt(:, :, 1:breaks(1)) = repmat(Theta(:, :, 1), [1, 1, breaks(1)]);
    for t = 2:m
        Theta_tt(:, :, breaks(t - 1) + 1 : breaks(t)) = repmat(Theta(:, :, t), [1, 1, breaks(t) - breaks(t - 1)]);
    end
    Theta_tt(:, :, breaks(m) + 1 : end) = repmat(Theta(:, :, m + 1), [1, 1, T - breaks(m)]);
else % if no break
    Y = mvnrnd(zeros(p,1), inv(Theta(:, :, 1)), T);
    Y_test = zeros(T, p, Nsim);
    for tt = 1 : Nsim
        Y_test(:, :, tt) = mvnrnd(zeros(p, 1), inv(Theta(:,:,1)), T);
    end
    Theta_tt = zeros(p, p, T);
    Theta_tt(:, :, 1:end) = repmat(Theta(:, :, 1), [1, 1, T]);
end
end