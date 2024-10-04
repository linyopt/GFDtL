function [AIC, BIC] = GetICs(Thetas, lossval)
%% Compute the AICs and BICs.
% AIC / BIC can be used to select hyperparameters.
% They are defined as
%       AIC = 2k - 2 ln(L),
%       BIC = k (ln(T) + ln(2pi)) - ln(L),
% where k is the number of estimated parameters; L is the
% maximized value of the likelihood function for the model,
% i.e., ln(L) is the negative log likelihood, the value of loss
% function; T is the sample size.
%
% - Usage:
%   [AIC, BIC] = GetICs(Thetas, lossvals);
%
% - Input:
%   @Thetas:    The obtained estimator.
%   @lossval:   The loss function value.
%
% - Output:
%   @AIC:       The AIC.
%   @BIC:       The BIC.

% Dimension and sample size.
[d, ~, T] = size(Thetas);

% The threshold for determining if an element is zero.
threshold_ = 1e-6;

% The indices of non-diagonal elements.
idx_diff = ~repmat(eye(d), [1, 1, T - 1]);

% The indices for the first matrix.
idx_1 = 1:d^2;

% The first-order difference of the input Theta.
diff_Theta = diff(Thetas, 1, 3);

% Compute K in the definition of AIC/BIC.
% It comtains two parts.
% 1. The number of non-zero differences of non-diagonal elements.
% 2. The number of non-zero elements of the non-diagonal elements of the
%    first matrix.
K = sum(abs(diff_Theta(idx_diff)) >= threshold_)...
    + sum(abs(Thetas(idx_1(~eye(d)))) >= threshold_, 'all');

% Compute AIC/BIC.
AIC = 2 * K + 2 * lossval;  
BIC = K * log(T) + lossval;
end