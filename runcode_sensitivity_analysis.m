%% Runcode for sensitivity analysis, i.e., Subsection 6.2 in the paper.
% Run the first section of this file with different `num_breaks` to get
% different columns of Figure 1 in the paper.
%
% It will take a very long time for each `num_breaks`.
%
% After finishing, you can only rerun the second/third section to obtain
% the figures (don't forget to change `num_breaks` accordingly).
%
% The second and third sections of this file are just two different ways to
% draw the figures.

%% Run experiment.
addpath(genpath(pwd))
clear
clc
% Settings
% The number of time slots.
T = 100;

% The dimension.
p = 10;

% The number of breaks.
num_breaks = 1; %0; %3;

% The sparsity.
% It can be a scalar or a vector. Normally one only need to give a positive
% scalar that is smaller than or equal to 1.
prob = 0.2;

% The number of numerical experiments.
num_exps = 1; 

% The candidates of lamb1s, lamb2s.
lamb1s = 0.1:0.01:1;
lamb2s = 10:1:200;

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
tic
[Y, Y_test, trueTheta, trueBreaks, Ests, breaksDTr, ThetaGs, cpGs, ...
 lossvals_DTr, lossvals_GFGL, BICs_DTr, BICs_GFGL, ...
 HDs_DTr, F1s_DTr, accs_DTr, errors_DTr, ...
 HDs_GFGL, F1s_GFGL, accs_GFGL, errors_GFGL, ...
 dresHFE, dparaHFE, dresBIC, dparaBIC, dresloss, dparaloss, ...
 gresHFE, gparaHFE, gresBIC, gparaBIC, gresloss, gparaloss, ...
 EstHFE, EstBIC, Estloss, ThGHFE, ThGBIC, ThGloss] = ...
er_exps(T, p, num_breaks, prob, num_exps, lamb1s, lamb2s, lamb3, epsilon, beta_, tau, Nsim, 0);
toc

filename = ['sensitivity-T', num2str(T), '-p', num2str(p), '-s', num2str(prob), '-m', num2str(num_breaks), '.mat'];
save(filename, 'Y', 'trueTheta', ...
    'lossvals_DTr', 'BICs_DTr', ...
    'HDs_DTr', 'F1s_DTr', 'accs_DTr', 'errors_DTr', ...
    'dresHFE', 'dparaHFE', 'dresBIC', 'dparaBIC', 'dresloss', 'dparaloss', ...
    'EstHFE', 'EstBIC', 'Estloss');

%% Draw "smooth" heatmap.

T = 100;
p = 10;
prob = 0.2;
num_breaks = 0;

filename = ['sensitivity-T', num2str(T), '-p', num2str(p), '-s', num2str(prob), ...
            '-m', num2str(num_breaks), '.mat'];
load(filename)

% The candidates of lamb1s, lamb2s.
lamb1s = 0.1:0.01:1;
lamb2s = 10:1:200;

figure
imagesc(lamb1s, lamb2s, HDs_DTr');
colorbar
colormap((sky(128)))
title("Hausdorff distances for D-trace loss (m^\ast=" + num_breaks + ")", 'FontSize', 14)
xlabel("\lambda_1", 'FontSize', 14)
ylabel("\lambda_2", 'FontSize', 14)
exportgraphics(gcf, "sen-ana-HDs-m" + num_breaks + ".png", 'Resolution', 300);

figure
imagesc(lamb1s, lamb2s, F1s_DTr');
colorbar
colormap(flipud(sky(128)))
title("F1 scores for D-trace loss (m^\ast=" + num_breaks + ")", 'FontSize', 14)
xlabel("\lambda_1", 'FontSize', 14)
ylabel("\lambda_2", 'FontSize', 14)
clim([0, 1])
exportgraphics(gcf, "sen-ana-F1s-m" + num_breaks + ".png", 'Resolution', 300);

figure
imagesc(lamb1s, lamb2s, accs_DTr');
colorbar
colormap(flipud(sky(128)))
title("Accuracies for D-trace loss (m^\ast=" + num_breaks + ")", 'FontSize', 14)
xlabel("\lambda_1", 'FontSize', 14)
ylabel("\lambda_2", 'FontSize', 14)
clim([0, 1])
exportgraphics(gcf, "sen-ana-accs-m" + num_breaks + ".png", 'Resolution', 300);

% figure
% imagesc(lamb1s, lamb2s, errors_DTr');
% colorbar
% colormap(flipud(flipud(sky(128))))
% title("Error for D-trace loss", 'FontSize', 14)
% xlabel("\lambda_1", 'FontSize', 14)
% ylabel("\lambda_2", 'FontSize', 14)

% Because BICs contain inf and the gap between maximum and minimum BICs is
% very large, directly drawing the heatmap provides nothing interesting.
% We need to process BICs.
% Make a copy.
BICs_DTr_ = BICs_DTr;
% Replace inf by the maximum (finite) value of BICs_DTr.
BICs_DTr_(isinf(BICs_DTr_)) = max(BICs_DTr(~isinf(BICs_DTr)));
% Minus the minimum value of BICs_DTr so that every value is nonnegative.
BICs_DTr_ = BICs_DTr_ - min(BICs_DTr_(:));
% Since all values are now nonnegative, we can apply log1p to all
% elements to get logarithmic scaling.
% Here we use log1p to handle 0.
BICs_DTr_ = log1p(BICs_DTr_);
figure
imagesc(lamb1s, lamb2s, BICs_DTr_');
colorbar
colormap(sky(128))
title("BICs for D-trace loss (m^\ast=" + num_breaks + ")", 'FontSize', 14)
xlabel("\lambda_1", 'FontSize', 14)
ylabel("\lambda_2", 'FontSize', 14)
exportgraphics(gcf, "sen-ana-BICs-m" + num_breaks + ".png", 'Resolution', 300);

% Because lossvals contain inf and the gap between maximum and minimum
% lossvals is very large, directly drawing the heatmap provides nothing
% interesting. 
% We need to process lossvals.
% Make a copy.
lossvals_DTr_ = lossvals_DTr;
% Replace inf by the maximum (finite) value of BICs_DTr.
lossvals_DTr_(isinf(lossvals_DTr_)) = max(lossvals_DTr(~isinf(lossvals_DTr)));
% Minus the minimum value of lossvals_DTr so that every value is nonnegative.
lossvals_DTr_ = lossvals_DTr_ - min(lossvals_DTr_(:));
% Since all values are now nonnegative, we can apply log1p to all
% elements to get logarithmic scaling.
% Here we use log1p to handle 0.
lossvals_DTr_ = log1p(lossvals_DTr_);
figure
imagesc(lamb1s, lamb2s, lossvals_DTr_');
colorbar
colormap(sky(128))
title("lossvals for D-trace loss (m^\ast=" + num_breaks + ")", 'FontSize', 14)
xlabel("\lambda_1", 'FontSize', 14)
ylabel("\lambda_2", 'FontSize', 14)
exportgraphics(gcf, "sen-ana-lossvals-m" + num_breaks + ".png", 'Resolution', 300);

%% Draw heatmap

T = 100;
p = 10;
prob = 0.2;
num_breaks = 0;

filename = ['sensitivity-T', num2str(T), '-p', num2str(p), '-s', num2str(prob), ...
            '-m', num2str(num_breaks), '.mat'];
load(filename)

% The candidates of lamb1s, lamb2s.
lamb1s = 0.1:0.01:1;
lamb2s = 10:1:200;

len1 = length(lamb1s);
len2 = length(lamb2s);

figure
Dh_HD = heatmap(lamb1s, lamb2s, HDs_DTr');
Dh_HD.GridVisible = 'off';
Dh_HD.Title = "Hausdorff distances for D-trace loss (m^\ast=" + num_breaks + ")";
Dh_HD.XLabel = "\lambda_1";
Dh_HD.YLabel = "\lambda_2";
Dh_HD.XDisplayLabels(mod(1:len1, 10) ~= 1) = {''};
Dh_HD.YDisplayLabels(mod(1:len2, 10) ~= 1) = {''};

figure
Dh_F1 = heatmap(lamb1s, lamb2s, F1s_DTr');
Dh_F1.GridVisible = 'off';
Dh_F1.Title = "F1 scores for D-trace loss (m^\ast=" + num_breaks + ")";
Dh_F1.XLabel = "\lambda_1";
Dh_F1.YLabel = "\lambda_2";
Dh_F1.XDisplayLabels(mod(1:len1, 10) ~= 1) = {''};
Dh_F1.YDisplayLabels(mod(1:len2, 10) ~= 1) = {''};
% Flip the colormap so that consistently lighter color is better.
dcmap = Dh_F1.Colormap;
Dh_F1.Colormap = flipud(dcmap);
clim([0, 1])

% figure
% Dh_error = heatmap(lamb1s, lamb2s, errors_DTr');
% Dh_error.GridVisible = 'off';
% Dh_error.Title = "Errors for D-trace loss (m^\ast=" + m + ")";
% Dh_error.XLabel = "\lambda_1";
% Dh_error.YLabel = "\lambda_2";
% Dh_error.XDisplayLabels(mod(1:len1, 10) ~= 1) = {''};
% Dh_error.YDisplayLabels(mod(1:len2, 10) ~= 1) = {''};

% Because BICs contain inf and the gap between maximum and minimum BICs is
% very large, directly drawing the heatmap provides nothing interesting.
% We need to process BICs.
% Make a copy.
BICs_DTr_ = BICs_DTr;
% Replace inf by the maximum (finite) value of BICs_DTr.
BICs_DTr_(isinf(BICs_DTr_)) = max(BICs_DTr(~isinf(BICs_DTr)));
% Minus the minimum value of BICs_DTr so that every value is nonnegative.
BICs_DTr_ = BICs_DTr_ - min(BICs_DTr_(:));
% Since all values are now nonnegative, we can apply log1p to all
% elements to get logarithmic scaling.
% Here we use log1p to handle 0.
BICs_DTr_ = log1p(BICs_DTr_);
figure
Dh_BICs = heatmap(lamb1s, lamb2s, BICs_DTr_');
Dh_BICs.GridVisible = 'off';
Dh_BICs.Title = "BICs for D-trace loss (m^\ast=" + num_breaks + ")";
Dh_BICs.XLabel = "\lambda_1";
Dh_BICs.YLabel = "\lambda_2";
Dh_BICs.XDisplayLabels(mod(1:len1, 10) ~= 1) = {''};
Dh_BICs.YDisplayLabels(mod(1:len2, 10) ~= 1) = {''};

% Because lossvals contain inf and the gap between maximum and minimum
% lossvals is very large, directly drawing the heatmap provides nothing
% interesting. 
% We need to process lossvals.
% Make a copy.
lossvals_DTr_ = lossvals_DTr;
% Replace inf by the maximum (finite) value of BICs_DTr.
lossvals_DTr_(isinf(lossvals_DTr_)) = max(lossvals_DTr(~isinf(lossvals_DTr)));
% Minus the minimum value of lossvals_DTr so that every value is nonnegative.
lossvals_DTr_ = lossvals_DTr_ - min(lossvals_DTr_(:));
% Since all values are now nonnegative, we can apply log1p to all
% elements to get logarithmic scaling.
% Here we use log1p to handle 0.
lossvals_DTr_ = log1p(lossvals_DTr_);
figure
Dh_loss = heatmap(lamb1s, lamb2s, lossvals_DTr_');
Dh_loss.GridVisible = 'off';
Dh_loss.Title = "lossvals for D-trace loss (m^\ast=" + num_breaks + ")";
Dh_loss.XLabel = "\lambda_1";
Dh_loss.YLabel = "\lambda_2";
Dh_loss.XDisplayLabels(mod(1:len1, 10) ~= 1) = {''};
Dh_loss.YDisplayLabels(mod(1:len2, 10) ~= 1) = {''};