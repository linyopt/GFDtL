%% T = 100 with 0 break
addpath(genpath(pwd))
clear
clc
rng(1,"twister");
% Settings
% The number of time slots.
T = 100;

% The dimension.
p = 10;

% The number of breaks.
num_breaks = 0;

% The number of numerical experiments.
num_exps = 100;

% The candidates of lamb1s, lamb2s.
lamb1s = 0.1:0.1:1;
lamb2s = 10:10:200;

% The lamb3.
lamb3 = 10;

% The epsilon.
epsilon = 0.01;

% The augmented parameter used in ADMM.
beta_ = 100;

% The factor used to rescale the Theta during DGP.
tau = 1;
% The number of simulated validation datasets.
Nsim = 10;

% Run the numerical experiments.
[Y, Y_test, trueTheta, trueBreaks, Ests, breaksDTr, ThetaGs, cpGs, ...
    lossvals_DTr, lossvals_GFGL, BICs_DTr, BICs_GFGL, ...
    HDs_DTr, F1s_DTr, accs_DTr, errors_DTr, ...
    HDs_GFGL, F1s_GFGL, accs_GFGL, errors_GFGL, ...
    dresHFE, dparaHFE, dresBIC, dparaBIC, dresloss, dparaloss, ...
    gresHFE, gparaHFE, gresBIC, gparaBIC, gresloss, gparaloss, ...
    EstHFE, EstBIC, Estloss, ThGHFE, ThGBIC, ThGloss] = ...
    banded_exps(T, p, num_breaks, num_exps, lamb1s, lamb2s, lamb3, epsilon, beta_, tau, Nsim);

% The following code computes the mean of performance metrics for all
% numerical experiments.
% Recall that it contains
%   HD, F1, error, accuracy, number of breaks.
mean(dresHFE)
mean(dresBIC)
mean(dresloss)
mean(gresHFE)
mean(gresBIC)
mean(gresloss)

filename = ['banded-T', num2str(T), '-p', num2str(p), '-m', num2str(num_breaks), '.mat'];
save(filename, 'Y', 'trueTheta', ...
    'lossvals_DTr', 'lossvals_GFGL', 'BICs_DTr', 'BICs_GFGL', ...
    'HDs_DTr', 'F1s_DTr', 'accs_DTr', 'errors_DTr', ...
    'HDs_GFGL', 'F1s_GFGL', 'accs_GFGL', 'errors_GFGL', ...
    'dresHFE', 'dparaHFE', 'dresBIC', 'dparaBIC', 'dresloss', 'dparaloss', ...
    'gresHFE', 'gparaHFE', 'gresBIC', 'gparaBIC', 'gresloss', 'gparaloss', ...
    'EstHFE', 'EstBIC', 'Estloss', 'ThGHFE', 'ThGBIC', 'ThGloss');

%% T = 150 with 1 break
addpath(genpath(pwd))
clear
clc
rng(2,"twister");
% Settings
% The number of time slots.
T = 150;

% The dimension.
p = 10;

% The number of breaks.
num_breaks = 1;

% The number of numerical experiments.
num_exps = 100;

% The candidates of lamb1s, lamb2s.
lamb1s = 0.1:0.1:1;
lamb2s = 10:10:200;

% The lamb3.
lamb3 = 10;

% The epsilon.
epsilon = 0.01;

% The augmented parameter used in ADMM.
beta_ = 100;

% The factor used to rescale the Theta during DGP.
tau = 1;
% The number of simulated validation datasets.
Nsim = 10;

% Run the numerical experiments.
[Y, Y_test, trueTheta, trueBreaks, Ests, breaksDTr, ThetaGs, cpGs, ...
    lossvals_DTr, lossvals_GFGL, BICs_DTr, BICs_GFGL, ...
    HDs_DTr, F1s_DTr, accs_DTr, errors_DTr, ...
    HDs_GFGL, F1s_GFGL, accs_GFGL, errors_GFGL, ...
    dresHFE, dparaHFE, dresBIC, dparaBIC, dresloss, dparaloss, ...
    gresHFE, gparaHFE, gresBIC, gparaBIC, gresloss, gparaloss, ...
    EstHFE, EstBIC, Estloss, ThGHFE, ThGBIC, ThGloss] = ...
    banded_exps(T, p, num_breaks, num_exps, lamb1s, lamb2s, lamb3, epsilon, beta_, tau, Nsim);


% The following code computes the mean of performance metrics for all
% numerical experiments.
% Recall that it contains
%   HD, F1, error, accuracy, number of breaks.
mean(dresHFE)
mean(dresBIC)
mean(dresloss)
mean(gresHFE)
mean(gresBIC)
mean(gresloss)

filename = ['banded-T', num2str(T), '-p', num2str(p), '-m', num2str(num_breaks), '.mat'];
save(filename, 'Y', 'trueTheta', ...
    'lossvals_DTr', 'lossvals_GFGL', 'BICs_DTr', 'BICs_GFGL', ...
    'HDs_DTr', 'F1s_DTr', 'accs_DTr', 'errors_DTr', ...
    'HDs_GFGL', 'F1s_GFGL', 'accs_GFGL', 'errors_GFGL', ...
    'dresHFE', 'dparaHFE', 'dresBIC', 'dparaBIC', 'dresloss', 'dparaloss', ...
    'gresHFE', 'gparaHFE', 'gresBIC', 'gparaBIC', 'gresloss', 'gparaloss', ...
    'EstHFE', 'EstBIC', 'Estloss', 'ThGHFE', 'ThGBIC', 'ThGloss');

%% T = 150 with 4 breaks
addpath(genpath(pwd))
clear
clc
rng(3,"twister");
% Settings
% The number of time slots.
T = 150;

% The dimension.
p = 10;

% The number of breaks.
num_breaks = 4;

% The number of numerical experiments.
num_exps = 100;

% The candidates of lamb1s, lamb2s.
lamb1s = 0.1:0.1:1;
lamb2s = 10:10:200;

% The lamb3.
lamb3 = 10;

% The epsilon.
epsilon = 0.01;

% The augmented parameter used in ADMM.
beta_ = 100;

% The factor used to rescale the Theta during DGP.
tau = 1;
% The number of simulated validation datasets.
Nsim = 10;

% Run the numerical experiments.
[Y, Y_test, trueTheta, trueBreaks, Ests, breaksDTr, ThetaGs, cpGs, ...
    lossvals_DTr, lossvals_GFGL, BICs_DTr, BICs_GFGL, ...
    HDs_DTr, F1s_DTr, accs_DTr, errors_DTr, ...
    HDs_GFGL, F1s_GFGL, accs_GFGL, errors_GFGL, ...
    dresHFE, dparaHFE, dresBIC, dparaBIC, dresloss, dparaloss, ...
    gresHFE, gparaHFE, gresBIC, gparaBIC, gresloss, gparaloss, ...
    EstHFE, EstBIC, Estloss, ThGHFE, ThGBIC, ThGloss] = ...
    banded_exps(T, p, num_breaks, num_exps, lamb1s, lamb2s, lamb3, epsilon, beta_, tau, Nsim);

% The following code computes the mean of performance metrics for all
% numerical experiments.
% Recall that it contains
%   HD, F1, error, accuracy, number of breaks.
mean(dresHFE)
mean(dresBIC)
mean(dresloss)
mean(gresHFE)
mean(gresBIC)
mean(gresloss)

filename = ['banded-T', num2str(T), '-p', num2str(p), '-m', num2str(num_breaks), '.mat'];
save(filename, 'Y', 'trueTheta', ...
    'lossvals_DTr', 'lossvals_GFGL', 'BICs_DTr', 'BICs_GFGL', ...
    'HDs_DTr', 'F1s_DTr', 'accs_DTr', 'errors_DTr', ...
    'HDs_GFGL', 'F1s_GFGL', 'accs_GFGL', 'errors_GFGL', ...
    'dresHFE', 'dparaHFE', 'dresBIC', 'dparaBIC', 'dresloss', 'dparaloss', ...
    'gresHFE', 'gparaHFE', 'gresBIC', 'gparaBIC', 'gresloss', 'gparaloss', ...
    'EstHFE', 'EstBIC', 'Estloss', 'ThGHFE', 'ThGBIC', 'ThGloss');
