# Break recovery in graphical networks with D-Trace loss

This repository contains the implementation of the **Group Fused D-trace LASSO (GFDtL)** estimator for identifying sparse precision matrices with breaks in time-series data.
The estimator is obtained by solving the following optimization problem:

$$
\min_{\Theta_t\succeq \epsilon I_p, 1 \leq t \leq T}\{\frac{1}{T}\sum_{t=1}^{T} \left[\text{tr}(\frac{1}{2}\Theta_t^2 X_tX^\top_t) - \text{tr}(\Theta_t)\right] + \lambda_{1T}\overset{T}{\underset{t=1}{\sum}}\underset{u\neq v}{\sum}|\Theta_{uv,t}|+\lambda_{2T}\overset{T}{\underset{t=2}{\sum}} \|\Theta_t - \Theta_{t-1}\|_F\}.
$$

where $\mathcal{X}\_T=(X_1,\ldots,X_T)$ is the sample, $\lambda_{1,T},\lambda_{2,T}$ are the tuning parameters.

## License

`GFDtL` is licensed under the MIT License (see [LICENSE](https://github.com/linyopt/GFDtL/blob/master/LICENSE)).

## Installation

1. Clone this repository.

   ```bash
   git clone https://github.com/linyopt/GFDtL.git
   cd GFDtL
   ```

2. Use Matlab with version *higher than R2023a* to run the files with name following the format `runcode_*.m`.

   *Note:* The implementation relies on functions [pagemtimes](https://www.mathworks.com/help/matlab/ref/pagemtimes.html) and [pageeig](https://www.mathworks.com/help/matlab/ref/pageeig.html) that are introduced in R2020a and R2023a, respectively, and [function argument validation](https://www.mathworks.com/help/matlab/ref/arguments.html) that is introduced in R2019b.
3. Normally, you should be able to run the code directly if your matlab version is higher than R2023a.

   If you find there are some dependence errors, please check requirements of three libraries in `GFGL/libs`.
   The folder `GFGL` is for the competing algorithm in [Regularized Estimation of Piecewise Constant Gaussian Graphical Models: The Group-Fused Graphical Lasso](https://www.tandfonline.com/doi/full/10.1080/10618600.2017.1302340), which depends on several toolboxes.
   Although we have tried our best to remove most of dependence, there may still be some omissions.

## Code structure

- The folder `@GFDtL` is the implementation of our algorithm, which is organized as a class.
  Interested readers can read `@GFDtL/GFDtL.m` for the declaration of properties/functions of the class, and read `@GFDtL/run.m` for the main content of our ADMM algorithm.
- The folder `simu` is the code for simulating data. Three main functions are `simu/DGP_Erdos_Renyi.m`, `simu/DGP_iid.m` and `simu/DGP_banded.m` for generating Erdos Renyi dataset, iid dataset and banded dataset, respectively.
  For more details, please read Section 6.1 in our paper.
- The folder `utils` contains some helper functions.
  The two files`hanning.m` and `perform_thresholding.m` are used to remove the toolbox dependence of the competing algorithm.
  The folder `boundedline` is used for plotting line with shaded error/confidence bounds, the credit goes for [kakearney/boundedline-pkg](https://github.com/kakearney/boundedline-pkg/).
  Please call  `help` for more details of these functions.
- There are six runcode files:
  1. `runcode_simu_er.m`: The runcode for experiments on simulated datasets simulated according to Setting (i) in Section 6.1 of the paper.
  2. `runcode_simu_banded.m`: The runcode for experiments on simulated datasets simulated according to Setting (ii) in Section 6.1 of the paper.
  3. `runcode_simu_iid.m`: The runcode for experiments on simulated datasets simulated according to Setting (iii) in Section 6.1 of the paper.
  4. `runcode_sensitivity_analysis.m` The runcode for empirical sensitivity analysis in Section 6.2 of the paper.
  5. `runcode_computational_complexity.m`: The runcode for empirical computational complexity analysis in Section 6.3 of the paper.
  6. `runcode_real.m`: The runcode for experiments on the real dataset `SP500.mat`.
     See Section 7 in our paper for more details about this dataset.

## Usage: Example of simulated dataset

```matlab
addpath(genpath(pwd));

T = 100; % The sample size.
p = 10; % The dimension.
m = 1; % The number of breaks.
% Generate Erdos Renyi data.
% Y is the training dataset, Y_test is the testing dataset, 
% Theta is the sequence of true precision matrices, breaks is the array of true breaks.
[Y, Y_test, Theta, breaks] = DGP_Erdos_Renyi(T, p, m);

% Declare GFDtL.
% Check `help GFDtL` for more details.
% Pass arguments by name, do NOT use positional arguments!
Est = GFDtL(Y=Y, lamb1=1, lamb2=50, lamb3=1, epsilon=0.01, beta_=21, disp_freq=1, maxiter=inf, tol=1e-3, tol_pcg=1e-2, tol_pcg_up=0.9);

% Run the algorithm.
Est.run

% Plot the difference of Theta.
% Check `help GFDtL.plot` for more details about the optional inputs.
Est.plot

% Estimate the breaks.
Est_breaks = Est.EstBreaks;
```

## Miscellaneous

Let `Est` be the GFDtL estimator.

### Possible unsolvability of the original problem

Our algorithm can detect the possible unsolvability of the original problem, which may occur when lamb2 is too small. If this happens, `Est.infeas` will be set to true.

### AIC / BIC

The **AIC/BIC** values for the estimator are stored in `Est.AIC` and `Est.BIC`. These are computed by calling `Est.GetICs`. If the algorithm terminates successfully, it will be called automatically, allowing you to access `Est.AIC` and `Est.BIC` directly. Otherwise, you can manually call `Est.GetICs` to compute these values.

### Built-in documents

See documents of functions by calling `help` in matlab for more details of the three DGP functions, `GFDtL`, and methods of `GFDtL`.