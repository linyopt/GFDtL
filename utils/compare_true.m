function [HD_DTr, F1_DTr, acc_DTr, error_DTr, HD_GFGL, F1_GFGL, acc_GFGL, error_GFGL] = ...
         compare_true(Est, ThetaGFGL, cpGs, trueTheta, trueBreaks)
%% Compare the obtained estimators and the true estimator, and compute the 
% metrics, including HD, error, F1.
%
% - Usage:
% [HD_DTr, F1_DTr, acc_DTr, error_DTr, HD_GFGL, F1_GFGL, acc_GFGL, error_GFGL] = ...
%          compare_true(Est, ThetaGFGL, cpGs, trueTheta, trueBreaks)
%
% - Input:
%   @Est:       The D-trace estimator.
%   @ThetaGFGL: The GFGL estimator.
%   @cpGs:      The break obtained by GFGL.
%   @trueTheta: The true Theta.
%   @trueBreaks:The true breaks.
%
% - Output:
%   @HD_DTr:    The HD for D-trace estimator.
%   @F1_DTr:    The F1 for D-trace estimator.
%   @acc_DTr:   The accuracy for D-trace estimator.
%   @error_DTr: The error for D-trace estimator.
%   @HD_GFGL:    The HD for GFGL estimator.
%   @F1_GFGL:    The F1 for GFGL estimator.
%   @acc_DTr:    The accuracy for GFGL estimator.
%   @error_GFGL: The error for GFGL estimator.

% Obtain the sample size T and the dimension p.
[~, p, T] = size(Est.Theta_k);

% Compute the number of all parameters for each matrix.
param = p^2;

% Estimate the breaks obtained by our model and algorithm, i.e., D-trace
% and ADMM.
breaksDTr = Est.EstBreaks();

% Compute HD, handle cases of one of the sets is empty.
% For non-break cases.
if isempty(trueBreaks)
    % If the break estimated by D-Trace is empty.
    if isempty(breaksDTr)
        % Set HD to be 0.
        HD_DTr = 0;
    else
        % Otherwise, add 0 to both sets.
        HD_DTr = max(breaksDTr); % = HausdorffDist([0], [0; breaksDTr]);
    end

    % If the break estimated by GFGL is empty.
    if isempty(cpGs)
        % Set HD to be 0.
        HD_GFGL = 0;
    else
        % Otherwise, add 0 to both sets.
        HD_GFGL = max(cpGs); % = HausdorffDist([0], [0; cpGs]);
    end
else
    % For normal cases.
    % If the break estimated by D-Trace is empty.
    if isempty(breaksDTr)
        % Add 0 to both sets.
        HD_DTr = max(trueBreaks); % = HausdorffDist([0; trueBreaks], [0; breaksDTr]);
    else
        % Otherwise, compute the HD.
        HD_DTr = HausdorffDist(trueBreaks, breaksDTr);
    end

    % If the break estimated by GFGL is empty.
    if isempty(cpGs)
        % Add 0 to both sets.
        HD_GFGL = max(trueBreaks); % = HausdorffDist([0; trueBreaks], [0; cpGs]);
    else
        % Otherwise, compute the HD.
        HD_GFGL = HausdorffDist(trueBreaks, cpGs);
    end
end

% Compute the averaged errors with respect to Theta.
error_DTr = sqrt(sum(vecnorm(reshape(Est.Theta_k - trueTheta, [], T)).^2) / param / T);
error_GFGL = sqrt(sum(vecnorm(reshape(ThetaGFGL - trueTheta, [], T)).^2) / param / T);

% Compute the F1 score.
F1_DTr = F1score(trueTheta, Est.Theta_k);
F1_GFGL = F1score(trueTheta, ThetaGFGL);

% Compute the accuracy.
acc_DTr = sum((trueTheta & abs(Est.Theta_k) >= 1e-5) | (~trueTheta & abs(Est.Theta_k) <= 1e-5), 'all') / ...
            (param * T);
acc_GFGL = sum((trueTheta & abs(ThetaGFGL) >= 1e-5) | (~trueTheta & abs(ThetaGFGL) <= 1e-5), 'all') / ...
            (param * T);
end