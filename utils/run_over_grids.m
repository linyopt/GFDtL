function [Ests, breaksDTr, ThetaGs, cpGs, ...
          lossvals_DTr, lossvals_GFGL, BICs_DTr, BICs_GFGL, ...
          HDs_DTr, F1s_DTr, accs_DTr, errors_DTr, ...
          HDs_GFGL, F1s_GFGL, accs_GFGL, errors_GFGL] = ...
run_over_grids(Y, Y_test, lamb1s, lamb2s, lamb3, epsilon, beta_, trueTheta, trueBreaks, gfgl)
%% Given candidates of lamb1s and lamb2s, obtain D-trace and GFGL estimators for
%  each pair of parameters.
%  Compute performance metrics and put them together.
%
% - Usage:
%   [Ests, breaksDTr, ThetaGs, cpGs, ...
%           lossvals_DTr, lossvals_GFGL, BICs_DTr, BICs_GFGL, ...
%           HDs_DTr, F1s_DTr, accs_DTr, errors_DTr, ...
%           HDs_GFGL, F1s_GFGL, accs_GFGL, errors_GFGL] = ...
% run_over_grids(Y, Y_test, lamb1s, lamb2s, lamb3, epsilon, beta_, trueTheta, trueBreaks, gfgl);
%
% - Input:
%   @Y:         The training dataset.
%   @Y_test:    The validation datasets.
%   @lamb1s:    The candidates of lamb1's.
%   @lamb2s:    The candidates of lamb2's.
%   @lamb3:     lamb3.
%   @epsilon:   epsilon.
%   @beta_:     The augmented parameter in ADMM.
%   @trueTheta: The true Theta.
%   @trueBreaks:The true breaks.
%   @gfgl:      The indicator for do we need to perform GFGL.
%
% - Output:
%   @Ests:      The D-trace estimators.
%   @breaksDTr: The breaks of D-trace estimators.
%   @ThetaGs:   The GFGL estimators.
%   @cpGs:      The breaks of GFGL estimators.
%   @lossvals_DTr:The loss values of D-trace estimators.
%   @lossvals_GFGL:The loss values of GFGL estimators.
%   @BICs_DTr:  The BICs of D-trace estimators.
%   @BICs_GFGL: The BICs of GFGL estimators.
%   @HDs_DTr:   The HDs of D-trace estimators.
%   @F1s_DTr:   The F1s of D-trace estimators.
%   @accs_DTr:  The accuracy of D-trace estimators.
%   @errors_DTr:The errors of D-trace estimators.
%   @HDs_GFGL:   The HDs of GFGL estimators.
%   @F1s_GFGL:   The F1s of GFGL estimators.
%   @accs_GFGL:  The accuracy of GFGL estimators.
%   @errors_GFGL:The errors of GFGL estimators.

if nargin == 9
    gfgl = 1;
end

tol = 1e-3;

% The number of candidate lamb1's.
len1 = length(lamb1s);
% The number of candidate lamb2's.
len2 = length(lamb2s);

% The cell of Estimators of DTrace loss.
Ests = cell(len1, len2);
% The call of breaks obtained by DTrace loss.
breaksDTr = cell(len1, len2);
% The loss values.
lossvals_DTr(len1, len2) = 0;
% The BICs.
BICs_DTr(len1, len2) = 0;
% The Hausdorff distances.
HDs_DTr(len1, len2) = 0;
% The F1 scores.
F1s_DTr(len1, len2) = 0;
% The accuracy.
accs_DTr(len1, len2) = 0;
% The errors.
errors_DTr(len1, len2) = 0;

% The cell of Theta's of GFGL.
ThetaGs = cell(len1, len2);
% The change points obtained by GFGL.
cpGs = cell(len1, len2);
% The loss values.
lossvals_GFGL(len1, len2) = 0;
% The BICs.
BICs_GFGL(len1, len2) = 0;
% The Hausdorff distances.
HDs_GFGL(len1, len2) = 0;
% The F1 scores.
F1s_GFGL(len1, len2) = 0;
% The accuracy.
accs_GFGL(len1, len2) = 0;
% The errors.
errors_GFGL(len1, len2) = 0;

parfor i = 1:len1
    for j = 1:len2
        % Our algorithm.
        Ests{i, j} = GFDtL(Y=Y, lamb1=lamb1s(i), lamb2=lamb2s(j), ...
                               lamb3=lamb3, epsilon=epsilon, beta_=beta_, ...
                               disp_freq=inf, maxiter=inf, tol=tol, ...
                               tol_pcg=1e-2, tol_pcg_up=0.9);
        Ests{i, j}.run;
        % Check infeasibility, if this choice of parameter causes the
        % infeasibility of the original D-trace problem, then let the loss
        % function value be inf so that this choice of parameters will be
        % never selected.
        if Ests{i, j}.infeas
            lossvals_DTr(i, j) = inf;
            BICs_DTr(i, j) = inf;
        else
            % Otherwise, compute the real loss function value.
            lossvals_DTr(i, j) = Ests{i, j}.lossfunc(Y_test);
            BICs_DTr(i, j) = Ests{i, j}.BIC;
        end
        breaksDTr{i, j} = Ests{i, j}.EstBreaks()';

        % GFGL.
        if gfgl
            [ThetaGs{i, j}, ~, cpGs{i, j}, ~, ~, ~] = GFGL(Y, lamb1s(i), ...
                                                           lamb2s(j), beta_);
            lossvals_GFGL(i, j) = lossfuncGau(ThetaGs{i, j}, Y_test);
            [~, BICs_GFGL(i, j)] = GetICs(ThetaGs{i, j}, lossfuncGau(ThetaGs{i, j}, Y));
    
            % Compare the obtained estimators with the ground truth and obtain
            % some metrics.
            [HDs_DTr(i, j), F1s_DTr(i, j), accs_DTr(i, j), errors_DTr(i, j), ...
             HDs_GFGL(i, j), F1s_GFGL(i, j), accs_GFGL(i, j), errors_GFGL(i, j)] = ...
             compare_true(Ests{i, j}, ThetaGs{i, j}, cpGs{i, j}, trueTheta, trueBreaks);
        else
            % Compare the obtained estimators with the ground truth and obtain
            % some metrics.
            [HDs_DTr(i, j), F1s_DTr(i, j), accs_DTr(i, j), errors_DTr(i, j), ...
             HDs_GFGL(i, j), F1s_GFGL(i, j), accs_GFGL(i, j), errors_GFGL(i, j)] = ...
             compare_true(Ests{i, j}, trueTheta, trueBreaks, trueTheta, trueBreaks);
        end
    end
end
end