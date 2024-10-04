function [Y, Y_test, trueTheta, trueBreaks, Ests, breaksDTr, ThetaGs, cpGs, ...
         lossvals_DTr, lossvals_GFGL, BICs_DTr, BICs_GFGL, ...
         HDs_DTr, F1s_DTr, accs_DTr, errors_DTr, ...
         HDs_GFGL, F1s_GFGL, accs_GFGL, errors_GFGL, ...
         dresHFE, dparaHFE, dresBIC, dparaBIC, dresloss, dparaloss, ...
         gresHFE, gparaHFE, gresBIC, gparaBIC, gresloss, gparaloss, ...
         EstHFE, EstBIC, Estloss, ThGHFE, ThGBIC, ThGloss] = ...
er_exps(T, p, num_breaks, prob, num_exps, lamb1s, lamb2s, lamb3, epsilon, beta_, tau, Nsim, gfgl)
%% Run experiments for Erdo-Renyi settings.
% Obtain estimators (both D-trace and GFGL) based on several simulated
% Erdo-Renyi datasets and compare the performance.
%
% - Usage:
%   [Y, Y_test, trueTheta, trueBreaks, Ests, breaksDTr, ThetaGs, cpGs, ...
%          lossvals_DTr, lossvals_GFGL, BICs_DTr, BICs_GFGL, ...
%          HDs_DTr, F1s_DTr, accs_DTr, errors_DTr, ...
%          HDs_GFGL, F1s_GFGL, accs_GFGL, errors_GFGL, ...
%          dresHFE, dparaHFE, dresBIC, dparaBIC, dresloss, dparaloss, ...
%          gresHFE, gparaHFE, gresBIC, gparaBIC, gresloss, gparaloss, ...
%          EstHFE, EstBIC, Estloss, ThGHFE, ThGBIC, ThGloss] = ...
% er_exps(T, p, num_breaks, prob, num_exps, lamb1s, lamb2s, lamb3, epsilon, beta_, tau, Nsim, gfgl);
%
% - Input:
%   @T:         The number of time slots; default 100.
%   @p:         The dimension; default 10.
%   @num_breaks:The number of breaks; default 0;
%   @prob:      The sparsity; default 0.2.
%   @num_exps:  The number of experiments; default 5.
%   @lamb1s:    The candidates of lamb1's; default 0.1:0.1:1.
%   @lamb2s:    The candidates of lamb2's; default 10:10:100.
%   @lamb3:     The lamb3; default 10.
%   @epislon:   The epsilon; default 0.01.
%   @beta_:     The augmented parameter used in ADMM; default 100.
%   @tau:       The scaling factor of Theta; default 1.
%   @Nsim:      The number of simulated validation datasets; default 10.
%   @gfgl:      The indicator for do we need to perform GFGL; default 1.
%
% - Output:
%   @Y:         The generated training datasets.
%   @Y_test:    The generated validation datasets.
%   @trueTheta: The true underlying Theta.
%   @trueBreaks:The true underlying breaks.
%   @Ests:      The obtained D-trace estimators.
%   @breaksDTr: The breaks obtained by D-trace estimators.
%   @ThetaGs:   The obtained GFGL estimators.
%   @cpGs:      The breaks obtained by GFGL estimators.
%   @lossvals_DTr:The values of loss function obtained by D-trace
%               estimators.
%   @lossvals_GFGL:The values of loss function obtained by GFGL estimators.
%   @BICs_DTr:  The BICs of D-trace estimators.
%   @BICs_GFGL: The BICs of GFGL estimators.
%   @HDs_DTr:   The Hausdorff distances of D-trace estimators.
%   @F1s_DTr:   The F1 scores of D-trace estimators.
%   @accs_DTr:  The accuracy of D-trace estimators.
%   @errors_DTr:The errors of D-trace estimators.
%   @HDs_GFGL:  The Hausdorff distance of GFGL estimators.
%   @F1s_GFGL:  The F1 scores of GFGL estimators.
%   @accs_GFGL: The accuracy of GFGL estimators.
%   @error_GFGL:The error s of GFGL estimators.
%   @dresHFE:   For each experiment, the best result for D-trace loss with
%               respect to the criterion (in short, HFE)
%                   minimum HD -> maximum F1 -> minimum error.
%   @dparaHFE:  For each experiment, the best tuning parameters for D-trace
%               loss with respect to the criterion (in short, HFE)
%                   minimum HD -> maximum F1 -> minimum error.
%   @dresBIC:   For each experiment, the best result for D-trace loss with
%               respect to the minimum BIC.
%   @dparaBIC:  For each experiment, the best tuning parameters for D-trace
%               loss with respect to the minimum BIC.
%   @dresloss:  For each experiment, the best result for D-trace loss with
%               respect to the minimum lossval.
%   @dparaloss: For each experiment, the best tuning parameters for D-trace
%               loss with respect to the minimum lossval.
%   @gresHFE:   For each experiment, the best result for GFGL with
%               respect to the criterion (in short, HFE)
%                   minimum HD -> maximum F1 -> minimum error.
%   @gparaHFE:  For each experiment, the best tuning parameters for GFGL
%               with respect to the criterion (in short, HFE)
%                   minimum HD -> maximum F1 -> minimum error.
%   @gresBIC:   For each experiment, the best result for GFGL with
%               respect to the minimum BIC.
%   @gparaBIC:  For each experiment, the best tuning parameters for GFGL
%               with respect to the minimum BIC.
%   @gresloss:  For each experiment, the best result for GFGL with
%               respect to the minimum lossval.
%   @gparaloss: For each experiment, the best tuning parameters for GFGL
%               with respect to the minimum lossval.
%   @EstHFE:    For each experiment, the best estimator for D-trace loss
%               with repsect to the criterion HFE.
%   @EstBIC:    For each experiment, the best estimator for D-trace loss
%               with repsect to the minimum BIC.
%   @Estloss:   For each experiment, the best estimator for D-trace loss
%               with repsect to the minimum lossval.
%   @ThGHFE:    For each experiment, the best estimator for GFGL
%               with repsect to the criterion HFE.
%   @ThGBIC:    For each experiment, the best estimator for GFGL
%               with repsect to the minimum BIC.
%   @ThGloss:   For each experiment, the best estimator for GFGL
%               with repsect to the minimum lossval.
%
% - Note:
%   1. prob can be a scale or a vector, normally you only need to pass a scale.
%   2. The returned best estimators are the matrices.

arguments (Input)
    T int32 = 100
    p int32 = 10
    num_breaks int32 = 0
    prob double = 0.2
    num_exps int32 = 5
    lamb1s double = 0.1:0.1:1
    lamb2s double = 10:10:100
    lamb3 double = 10
    epsilon double = 0.01
    beta_ double = 100
    tau double = 1
    Nsim int32 = 10
    gfgl int32 = 1
end

% The number of candidates of lamb1, lamb2.
len1 = length(lamb1s);
len2 = length(lamb2s);

% The data for training.
Y(T, p, num_exps) = 0;

% The data for cross validation.
% There will be NSim datasets that are out of the training set but generated
% using the same underlying groundtruth.
% The the mean of values of loss function over the NSim datasets will be
% one of the criterions for us to select the lamb1, lamb2.
Y_test(T, p, Nsim, num_exps) = 0;

% The arrays for storage of true Theta's.
trueTheta(p, p, T, num_exps) = 0;
% The arrars for storage of true breaks.
if num_breaks ~= 0
    trueBreaks(num_breaks, num_exps) = 0;
end

% The cell of all Estimators of DTrace loss.
Ests = cell(num_exps, 1);
% The cell of all breaks estimated by DTrace loss.
breaksDTr = cell(num_exps, 1);
% The cell of all estimated Theta of GFGL.
ThetaGs = cell(num_exps, 1);
% The cell of all breaks estimated by GFGL.
cpGs = cell(num_exps, 1);
% All loss values.
lossvals_DTr = zeros(len1, len2, num_exps);
lossvals_GFGL = zeros(len1, len2, num_exps);
% All BICs.
BICs_DTr = zeros(len1, len2, num_exps);
BICs_GFGL = zeros(len1, len2, num_exps);
% HDs, F1s, accs, errors.
HDs_DTr = zeros(len1, len2, num_exps);
F1s_DTr = zeros(len1, len2, num_exps);
accs_DTr = zeros(len1, len2, num_exps);
errors_DTr = zeros(len1, len2, num_exps);
HDs_GFGL = zeros(len1, len2, num_exps);
F1s_GFGL = zeros(len1, len2, num_exps);
accs_GFGL = zeros(len1, len2, num_exps);
errors_GFGL = zeros(len1, len2, num_exps);

% The optimal performance results selected via different criterions.
% 1. Minimum HD -> maximum F1 -> minimum error for D-trace loss.
dresHFE = zeros(num_exps, 5);
dparaHFE = zeros(num_exps, 2);
% 2. Minimum HD -> maximum F1 -> minimum error for GFGL.
gresHFE = zeros(num_exps, 5);
gparaHFE = zeros(num_exps, 2);
% 3. Minimum BIC for D-trace loss.
dresBIC = zeros(num_exps, 5);
dparaBIC = zeros(num_exps, 2);
% 4. Minimum BIC for GFGL.
gresBIC = zeros(num_exps, 5);
gparaBIC = zeros(num_exps, 2);
% 5. Minimum lossval for D-trace loss.
dresloss = zeros(num_exps, 5);
dparaloss = zeros(num_exps, 2);
% 6. Minimum lossval for GFGL.
gresloss = zeros(num_exps, 5);
gparaloss = zeros(num_exps, 2);

% The optimal estimators selected via different criterions.
% 1. Minimum HD -> maximum F1 -> minimum error for D-trace loss.
EstHFE = cell(num_exps, 1);
% 2. Minimum HD -> maximum F1 -> minimum error for GFGL.
ThGHFE = cell(num_exps, 1);
% 3. Minimum BIC for D-trace loss.
EstBIC = cell(num_exps, 1);
% 4. Minimum BIC for GFGL.
ThGBIC = cell(num_exps, 1);
% 5. Minimum lossval for D-trace loss.
Estloss = cell(num_exps, 1);
% 6. Minimum lossval for GFGL.
ThGloss = cell(num_exps, 1);

for it = 1:num_exps
    fprintf("Experiment %g \n", it);
    % Call the simulating function to get the training dataset and
    % validating datasets.
    [Y(:, :, it), Y_test(:, :, :, it), trueTheta(:, :, :, it), ...
            trueBreaks(:, it)] = DGP_Erdos_Renyi(T, p, num_breaks, prob, tau, Nsim);

    % Run algorithms (D-trace and GFGL) over grids of lamb1/lamb2.
    [Ests{it}, breaksDTr{it}, ThetaGs{it}, cpGs{it}, lossvals_DTr(:, :, it), ...
     lossvals_GFGL(:, :, it), BICs_DTr(:, :, it), BICs_GFGL(:, :, it), ...
     HDs_DTr(:, :, it), F1s_DTr(:, :, it), accs_DTr(:, :, it), errors_DTr(:, :, it), ...
     HDs_GFGL(:, :, it), F1s_GFGL(:, :, it), accs_GFGL(:, :, it), errors_GFGL(:, :, it)] = ...
    run_over_grids(Y(:, :, it), Y_test(:, :, :, it), lamb1s, lamb2s, lamb3, ...
                epsilon, beta_, trueTheta(:, :, :, it), trueBreaks(:, it), gfgl);

    % 1. Minimum HD -> maximum F1 -> minimum error for D-trace loss.
    [dresHFE(it, :), dparaHFE(it, :), EstHFE{it}] = findHFE(HDs_DTr(:, :, it), ...
        F1s_DTr(:, :, it), errors_DTr(:, :, it), accs_DTr(:, :, it), breaksDTr{it}, ...
        lamb1s, lamb2s, Ests{it});
    EstHFE{it} = EstHFE{it}.Theta_k;

    % 2. Minimum HD -> maximum F1 -> minimum error for GFGL.
    [gresHFE(it, :), gparaHFE(it, :), ThGHFE{it}] = findHFE(HDs_GFGL(:, :, it), ...
        F1s_GFGL(:, :, it), errors_GFGL(:, :, it), accs_GFGL(:, :, it), cpGs{it}, ...
        lamb1s, lamb2s, ThetaGs{it});

    % 3. Minimum BIC for D-trace loss.
    [dresBIC(it, :), dparaBIC(it, :), EstBIC{it}] = findMin(BICs_DTr(:, :, it), ...
        HDs_DTr(:, :, it), F1s_DTr(:, :, it), errors_DTr(:, :, it), accs_DTr(:, :, it), ...
        breaksDTr{it}, lamb1s, lamb2s, Ests{it});
    EstBIC{it} = EstBIC{it}.Theta_k;

    % 4. Minimum BIC for GFGL.
    [gresBIC(it, :), gparaBIC(it, :), ThGBIC{it}] = findMin(BICs_GFGL(:, :, it), ...
        HDs_GFGL(:, :, it), F1s_GFGL(:, :, it), errors_GFGL(:, :, it), accs_GFGL(:, :, it), ...
        cpGs{it}, lamb1s, lamb2s, ThetaGs{it});

    % 5. Minimum lossval for D-trace loss.
    [dresloss(it, :), dparaloss(it, :), Estloss{it}] = findMin(lossvals_DTr(:, :, it), ...
        HDs_DTr(:, :, it), F1s_DTr(:, :, it), errors_DTr(:, :, it), accs_DTr(:, :, it), ...
        breaksDTr{it}, lamb1s, lamb2s, Ests{it});
    Estloss{it} = Estloss{it}.Theta_k;

    % 6. Minimum lossval for GFGL.
    [gresloss(it, :), gparaloss(it, :), ThGloss{it}] = findMin(lossvals_GFGL(:, :, it), ...
        HDs_GFGL(:, :, it), F1s_GFGL(:, :, it), errors_GFGL(:, :, it), accs_GFGL(:, :, it), ...
        cpGs{it}, lamb1s, lamb2s, ThetaGs{it});
end
end

function [resHFE, paraHFE, Est] = findHFE(HDs, F1s, errors, accuracy, breaks, lamb1s, lamb2s, Ests)
%% Find minimum HD -> maximum F1 -> minimum error.

% The result corresponding to this criterion, it will contain
%   HD, F1, error, accuracy, number of breaks.
resHFE = zeros(1, 5);

% Step 1: Find the minimum value in HDs and its indices.
resHFE(1) = min(HDs(:));
indHD = find(HDs == resHFE(1));

% Step 2: From the indHD, find the maximum value in F1s.
[resHFE(2), indF1InHD] = max(F1s(indHD));
indF1 = indHD(indF1InHD);

% Step 3: From the indF1, find the minimum value in errors.
[resHFE(3), indEInF1] = min(errors(indF1));
indE1 = indF1(indEInF1);

% Return the first one to avoid there are multiple results.
% This will be the obtained index.
idx = indE1(1);
% Accuracy.
resHFE(4) = accuracy(idx);
% Number of breaks.
resHFE(5) = length(breaks{idx});

% The corresponding tuning parameters.
[r, c] = ind2sub(size(HDs), idx);
paraHFE = [lamb1s(r), lamb2s(c)];

% The corresponding estimators.
Est = Ests{idx};
end

function [res, para, Est] = findMin(vals, HDs, F1s, errors, accuracy, breaks, lamb1s, lamb2s, Ests)
%% Find minimum BIC/lossval.

% Find the minimum value.
idx = find(vals == min(vals(:)));
% Return the first one to avoid there are multiple results.
% This will be the obtained index.
idx = idx(1);

% The result corresponding the this criterion, it will contain
%   HD, F1, error, accuracy, number of breaks.
res = [HDs(idx), F1s(idx), errors(idx), accuracy(idx), length(breaks{idx})];

% The corresponding tuning parameters.
[r, c] = ind2sub(size(HDs), idx);
para = [lamb1s(r), lamb2s(c)];

% The corresponding estimators.
Est = Ests{idx};
end