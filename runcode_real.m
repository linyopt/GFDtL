%% Real data 1.
addpath(genpath(pwd));
clear
clc
close all

load 'SP500.mat'
% % The columns contain (from left to right) the following stock indices:
% % Berkshire, JPMorgan, Exxon, General Electric, Apple, Lockheed Martin,
% % Procter&Gamble, Walmart, Amazon, Goldman Sachs, American Airlines
% % Alphabet, UPS, Verizon, Boeing, Chevron, Ford, Pfizer, Equity Residential
% % Jacobs Engineering Group

data = SP500.data; 
dates = SP500.dates; 
Y = 100*(log(data(2:end,:))-log(data(1:end-1,:)));

% Parameters and settings.
[T,d] = size(Y);

lamb1s = 0.1:0.1:2;
lamb2s = 10:10:200;
lamb3 = 50;
epsilon = 0;
beta_ = 100;
len1 = length(lamb1s);
len2 = length(lamb2s);

BIC_DTr = zeros(len1, len2);
BIC_GFGL = zeros(len1, len2);
breaks_DTr = cell(len1, len2);
breaks_GFGL = cell(len1, len2);

Ests = cell(len1, len2);
ThetaGs = cell(len1, len2);

parfor i = 1:len1
    for j = 1:len2
        Ests{i, j} = GFDtL(Y=Y, lamb1=lamb1s(i), lamb2=lamb2s(j), ...
                               lamb3=lamb3, epsilon=epsilon, beta_=beta_, ...
                               disp_freq=inf, maxiter=1e5, tol=1e-3, ...
                               tol_pcg=1e-2, tol_pcg_up=0.9);
        Ests{i, j}.run;
        % Compute BIC.
        BIC_DTr(i, j) = Ests{i, j}.BIC;
        % Compute break.
        breaks_DTr{i, j} = Ests{i, j}.EstBreaks();

        % Run GFGL.
        [ThetaGs{i, j}, ~, breaks_GFGL{i, j}, ~, ~, ~] = GFGL(Y, lamb1s(i), ...
                                                              lamb2s(j), beta_);

        % Compute BIC.
        lossG = lossfuncGau(ThetaGs{i, j}, Y);
        [~, BIC_GFGL(i, j)] = GetICs(ThetaGs{i, j}, lossG);
    end
end

% Post-process to obtain the final estimator by selecting the minimum BIC
% or the minimal lossval.
% Note that there may be more than one combination of parameters that can
% obtain the minimum BIC.
% Here we only pick the first combination of parameters.

% Obtain the minimum BIC.
[min_BIC_D, ind_BIC_D] = min(BIC_DTr, [], 'all');
min_BIC_Est = Ests{ind_BIC_D};
min_BIC_breaks_D = breaks_DTr{ind_BIC_D};

[min_BIC_G, ind_BIC_G] = min(BIC_GFGL, [], 'all');
min_BIC_ThetaG = ThetaGs{ind_BIC_G};
min_BIC_breaks_G = breaks_DTr{ind_BIC_G};

% Save BICs and breaks for all choices of parameters.
save("real_SP500.mat", "BIC_DTr", "BIC_GFGL", "breaks_DTr", "breaks_GFGL", ...
    "min_BIC_Est", "min_BIC_ThetaG", "min_BIC_breaks_G");

% Estimators.
Theta_d = min_BIC_Est.Theta_k;
Theta_g = min_BIC_ThetaG;

w_d = zeros(d,T); e_d = zeros(T,1);
w_g = zeros(d,T); e_g = zeros(T,1);
e_equi = zeros(d,1);

for t = 2:T
    w_d(:,t)= GMVP(inv(Theta_d(:,:,t-1))); w_g(:,t)= GMVP(inv(Theta_g(:,:,t-1)));
    e_d(t) = w_d(:,t)'*Y(t,:)'; e_g(t) = w_g(:,t)'*Y(t,:)'; 
    e_equi(t) = sum(Y(t,:))/d;
end
e_d = e_d(2:end); e_g = e_g(2:end); e_equi = e_equi(2:end); 

% out-of-sample GMVP returns
e_gmvp = [e_d e_g e_equi];

% out-of-sample average portfolio returns, standard deviations and
% information ratios
Results = [252*mean(e_gmvp);sqrt(252)*std(e_gmvp);(252*mean(e_gmvp))./(sqrt(252)*std(e_gmvp))];
