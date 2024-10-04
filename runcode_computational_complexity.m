%% Runcode for empirical computational comlexity analysis, i.e., Section 6.3 in the paper.
% Directly run the WHOLE file to get the four figures used in the paper.
%
% It may take a long time.
%
% After finishing, you can only rerun the final section of this file to get
% the figures.

%% Influence of sample size T
addpath(genpath(pwd));
clear
clc
close all

% Fix dimension p = 10; number of breaks = 1; prob = 0.2.
% We consider the Erdos Renyi DGP.
p = 10;
num_breaks = 1;
prob = 0.2;

% Algorithm parameters.
% For different T's, we fixed the same lamb1/lamb2, lmab3, epsilon, beta_.
lamb1 = 1;
lamb2 = 40;
lamb3 = 10;
epsilon = 0;
beta_ = 21;

% The candidates of T's.
Ts = [20, 50, 100, 150, 200, 250, 300];

% Time
times_T = zeros(length(Ts), 5);

for i = 1:length(Ts)
    for j = 1:5
        % Generate datasets.
        [Y, ~, ~, ~] = DGP_Erdos_Renyi(Ts(i), p, num_breaks, prob, 1, 0);
        % Run our algorithm and record time.
        tic
        Est = GFDtL(Y=Y, lamb1=lamb1, lamb2=lamb2, lamb3=lamb3, epsilon=epsilon, ...
                        beta_=beta_, disp_freq=inf, maxiter=inf, tol=1e-3, ...
                        tol_pcg=1e-2, tol_pcg_up=0.9);
        Est.run;
        times_T(i, j) = toc;
    end
end

%% Influence of dimension p
addpath(genpath(pwd))

% Fix sample size T = 100; number of breaks = 1; prob = 0.2.
% We consider the Erdos Renyi DGP.
T = 100;
num_breaks = 1;
prob = 0.2;

% Algorithm parameters.
% For different p's, we fixed the same lamb1/lamb2, lamb3, epsilon, beta_.
lamb1 = 1;
lamb2 = 40;
lamb3 = 10;
epsilon = 0;
beta_ = 21;

% The candidates of p's.
ps = [3, 6, 10, 30, 40, 60];

% Time
times_p = zeros(length(ps), 5);

for i = 1:length(ps)
    for j = 1:5
        % Generate datasets.
        [Y, ~, ~, ~] = DGP_Erdos_Renyi(T, ps(i), num_breaks, prob, 1, 0);
        % Run our algorithm and record time.
        tic
        Est = GFDtL(Y=Y, lamb1=lamb1, lamb2=lamb2, lamb3=lamb3, epsilon=epsilon, ...
                        beta_=beta_, disp_freq=inf, maxiter=inf, tol=1e-3, ...
                        tol_pcg=1e-2, tol_pcg_up=0.9);
        Est.run;
        times_p(i, j) = toc;
    end
end

%% Influence of the number of breaks
addpath(genpath(pwd));

% Fix sample size T = 100; dimension p = 10; prob = 0.2.
% We consider the Erdos Renyi DGP.
T = 100;
p = 10;
prob = 0.2;

% Algorithm parameters.
% For different numbers of breaks, we fixed the same lamb1/lamb2, lamb3, epsilon, beta_.
lamb1 = 1;
lamb2 = 40;
lamb3 = 10;
epsilon = 0;
beta_ = 21;

% The candidates of numbers of breaks.
ms = [0, 1, 3, 5, 7, 15, 30, 50];

% Time
times_m = zeros(length(ms), 5);

for i = 1:length(ms)
    for j = 1:5
        % Generate datasets.
        [Y, ~, ~, ~] = DGP_Erdos_Renyi(T, p, ms(i), prob, 1, 0, 'constant');
        % Run our algorithm and record time.
        tic
        Est = GFDtL(Y=Y, lamb1=lamb1, lamb2=lamb2, lamb3=lamb3, epsilon=epsilon, ...
                        beta_=beta_, disp_freq=inf, maxiter=inf, tol=1e-3, ...
                        tol_pcg=1e-2, tol_pcg_up=0.9);
        Est.run;
        times_m(i, j) = toc;
    end
end

%% Influence of lamb1/lamb2
addpath(genpath(pwd));

% Fix sample size T = 100; dimension p = 10; number of breaks = 1; prob = 0.2.
% We consider the Erdos Renyi DGP.
T = 100;
p = 10;
num_breaks = 1;
prob = 0.2;

% Generate datasets.
[Y, ~, ~, ~] = DGP_Erdos_Renyi(T, p, num_breaks, prob, 1, 0);

% Algorithm parameters.
% For different numbers of breaks, we fixed the same lamb3, epsilon, beta_.
lamb3 = 10;
epsilon = 0;
beta_ = 21;

% The candidates lamb1, lamb2.
lamb1s = 0.1:0.1:1;
lamb2s = 10:10:200;

len1 = length(lamb1s);
len2 = length(lamb2s);

% Time
times_ls = zeros(len1, len2, 5);

for i = 1:len1
    for j = 1:len2
        for k = 1:5
            % Run our algorithm and record time.
            tic
            Est = GFDtL(Y=Y, lamb1=lamb1s(i), lamb2=lamb2s(j), lamb3=lamb3, epsilon=epsilon, ...
                            beta_=beta_, disp_freq=inf, maxiter=inf, tol=1e-3, ...
                            tol_pcg=1e-2, tol_pcg_up=0.9);
            Est.run;
            times_ls(i, j, k) = toc;
        end
    end
end

times_ls = mean(times_ls, 3);

%% Influence of lamb3.

addpath(genpath(pwd));

% Fix sample size T = 100; dimension p = 10; number of breaks = 1; prob = 0.2.
% We consider the Erdos Renyi DGP.
T = 100;
p = 10;
num_breaks = 1;
prob = 0.2;

% Algorithm parameters.
% For different numbers of breaks, we fixed the same lamb1/lamb2, epsilon, beta_.
lamb1 = 1;
lamb2 = 40;
epsilon = 0;
beta_ = 21;

% The candidates lamb3.
lamb3s = [0.5, 1, 5, 10, 20, 30, 100, 150, 200];

% Time
times_l3s = zeros(length(lamb3s), 5);

for j = 1:5
    % Generate datasets.
    [Y, ~, ~, ~] = DGP_Erdos_Renyi(T, p, num_breaks, prob, 1, 0);
    for i = 1:length(lamb3s)
        % Run our algorithm and record time.
        tic
        Est = GFDtL(Y=Y, lamb1=lamb1, lamb2=lamb2, lamb3=lamb3s(i), epsilon=epsilon, ...
                        beta_=beta_, disp_freq=inf, maxiter=inf, tol=1e-3, ...
                        tol_pcg=1e-2, tol_pcg_up=0.9);
        Est.run;
        times_l3s(i, j) = toc;
    end
end

% Save.
save("times.mat", ...
     "Ts", "ps", "ms", "lamb1s", "lamb2s", "lamb3s", ...
     "times_T", "times_p", "times_m", "times_ls", "times_l3s");

%% Load data and draw figures.

load("times.mat");

% T.
figure();
mean_T = mean(times_T, 2);
std_T = std(times_T, 0, 2);
boundedline(Ts, mean_T, std_T, 'x-', 'alpha');
xlabel("Sample size T", 'FontSize', 14);
ylabel("Computation time (s)", 'FontSize', 14);
title("Computation time v.s. T (p=10, m^*=1, \lambda_1=1, \lambda_2=40, \lambda_3=10)", "FontSize", 14);
exportgraphics(gcf, 'comp-time-T.png', 'Resolution', 300);

% p.
figure();
mean_p = mean(times_p, 2);
std_p = std(times_p, 0, 2);
boundedline(ps, mean_p, std_p, 'x-', 'alpha');
xlabel("Dimension p", 'FontSize', 14);
ylabel("Computation time(s)", 'FontSize', 14)
title("Computation time v.s. p (T=100, m^*=1, \lambda_1=1, \lambda_2=40, \lambda_3=10)", 'FontSize', 14);
exportgraphics(gcf, 'comp-time-p.png', 'Resolution', 300);

% number of breaks.
figure();
mean_m = mean(times_m, 2);
std_m = std(times_m, 0, 2);
boundedline(ms, mean_m, std_m, 'x-', 'alpha');
xlabel("Number of true breaks", 'FontSize', 14);
ylabel("Computation time (s)", 'FontSize', 14)
title("Computation time v.s. m^* (T=100, p=10, \lambda_1=1, \lambda_2=40, \lambda_3=10)", 'FontSize', 14);
exportgraphics(gcf, 'comp-time-m.png', 'Resolution', 300);

% lamb1s / lamb2s.
figure();
h_ls = heatmap(lamb1s, lamb2s, times_ls');
h_ls.GridVisible = 'off';
h_ls.Title = "Computation time v.s. lamb1/lamb2 (T=100, p=10, m^*=1, \lambda_3=10)";
h_ls.XLabel = "\lambda_1";
h_ls.YLabel = "\lambda_2";
exportgraphics(gcf, 'comp-time-ls.png', 'Resolution', 300);

% lamb3s.
figure();
mean_l3 = mean(times_l3s, 2);
std_l3 = std(times_l3s, 0, 2);
boundedline(lamb3s, mean_l3, std_l3, 'x-', 'alpha');
xlabel("\lambda_3", 'FontSize', 14);
ylabel("Computation time (s)", 'FontSize', 14)
title("Computation time v.s. \lambda_3 (T=100, p=10, m^*=1, \lambda_1=1, \lambda_2=40)", 'FontSize', 14);
exportgraphics(gcf, 'comp-time-l3.png', 'Resolution', 300);