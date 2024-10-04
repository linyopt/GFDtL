classdef GFDtL < handle
    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% The data-related stuff.                                     %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The sample size.
        T double {mustBePositive, mustBeInteger}
        
        % The dimension of the vector in each time slot.
        d double {mustBePositive, mustBeInteger}
        
        % The time series data of size (d, T).
        Y double {mustBeNumeric}
        % The matrixized data of size (1, d, T).
        Ym double {mustBeNumeric}
        % The transport of matrixized data of size (d, 1, T).
        Ymt double {mustBeNumeric}
        
        % The covariance matrices obtained by Y. The size of Sigma is (d, d, T).
        Sigma double {mustBeNumeric}

        % The sqrt of Sigma.
        SqrtSigma double {mustBeNumeric}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% The hyperparameters and settings.                           %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The hyperparameter lamb1 for the lasso regularizer.
        % Default 0.1.
        lamb1 (1, 1) double {mustBeNonnegative} = 0.1
        
        % The hyperparameter lamb2 for the fused lasso regularizer.
        % Default 100.
        lamb2 (1, 1) double {mustBePositive} = 100
        
        % The hyperparameter lamb3 used in the modified fused lasso regularizer.
        % It should be >= 0.5 to ensure the convexity.
        % Default 0.5.
        lamb3 (1, 1) double {mustBeGreaterThanOrEqual(lamb3, 0.5)} = 0.5
        
        % The hyperparameter epsilon used for the positive definite constraint
        % in the D-trace loss. It should be positive.
        % Default 0.01
        epsilon (1, 1) double {mustBeNonnegative} = 0.01
        
        % The augmented Lagrangian hyperparameter beta_ to balance the primal
        % and dual. It should be positive.
        % Default 1.
        beta_ (1, 1) double {mustBePositive} = 1
        
        % The tolerance to stop the algorithm.
        % Default 1e-2.
        tol (1, 1) double {mustBePositive, mustBeLessThan(tol, 1)} = 1e-2
        
        % The initial tolerance for pcg to update Theta in each step.
        % Default 1e-3.
        tol_pcg (1, 1) double {mustBePositive, mustBeLessThan(tol_pcg, 1)} = 1e-3
        
        % The parameter to update tol_pcg.
        % For each 10 iterations, tol_pcg will decrease by multiplying tol_pcg_up.
        % So, tol_pcg_up has to be positive and strictly less than 1.
        % Default 0.7.
        tol_pcg_up (1, 1) double {mustBePositive, mustBeLessThanOrEqual(tol_pcg_up, 1)} = 0.7
        
        % Maximal iterations.
        % Default 5000.
        maxiter (1, 1) double {mustBeNumeric, mustBePositive} = 5000
        
        % Display frequency.
        % Default 1.
        disp_freq (1, 1) double {mustBePositive} = 1

        % Show figure or not.
        % Default false.
        showfig (1, 1) = false

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% The followings are declarations of primal / dual variables. %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Primal variables.
        % Theta
        Theta_k double
        % V
        V_k double
        % Upsilon
        Upsilon_k double
        % D_k
        D_k double
        
        % Dual variables.
        % Y
        Y_k double
        % Z
        Z_k double
        % A
        A_k double
        % Auxiliary variable zeta_.
        zeta_ double

        % Norm of first-order difference of Theta.
        norm_diff double

        % Helper variable indicating indices of diagonal elements.
        idx

        % Intermediate variables for acceleration.
        first_col
        third_col
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% The followings are lossfunc / primval / dualval.            %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The loss function value.
        lossval (1, 1) double = inf
        
        % The primal objective value.
        primval (1, 1) double = inf
        
        % The dual objective value.
        dualval (1, 1) double = inf

        % The AIC value.
        AIC (1, 1) double = inf

        % The BIC value.
        BIC (1, 1) double = inf

        % breaks
        breaks double

        % Infeasibility
        infeas logical = false
    end
    methods
        function obj = GFDtL(opts)
            %% The estimator based on D-trace loss.
            % - Usage:
            %   Est = GFDtL(Y=Y, lamb1=lamb1, lamb2=lamb2, lamb3=lamb3, 
            %                   beta_=beta_, epsilon=epsilon, tol=tol, 
            %                   tol_pcg=tol_pcg, tol_pcg_up=tol_pcg_up, 
            %                   maxiter=maxiter, disp_freq=disp_freq, 
            %                   showfig=showfig);
            % - Note:
            %   1. Please declare arguments in an explicit way like "Y=Y", 
            %      do not use positional arguments.
            %   2. The order of input arguments does not matter, i.e., 
            %      the following two ways are equivalent:
            %      >> Est = GFDtL(Y=Y, lamb1=lamb1);
            %      >> Est = GFDtL(lamb1=lamb1, Y=Y);
            %   3. Just skip the arguments that you want it to use the default
            %      value. For example, skipping "lamb2" will let "lamb2" be
            %      1e-1 by default.
            %   4. Normally, you only need to provide
            %   Y, lamb1, lamb2, lamb3, epsilon, tol, beta_, maxiter, disp_freq.
            %
            % - Input:
            %   @Y:         Data, a Txp matrix.
            %   @lamb1:     lamb1, positive scalar, default 1e-3.
            %   @lamb2:     lamb2, positive scalar, default 1e-1.
            %   @lamb3:     lamb3, positive scalar >= 0.5, default 0.5.
            %   @beta_:     Argumented Lagrangian parameter, positive scalar,
            %               default 1.
            %   @epsilon:   Parameter for PSD constraint, positive scalar, 
            %               default 0.01.
            %   @tol:       Tolerance to stop the algorithm, positive scalar < 1, 
            %               default 1e-2.
            %   @tol_pcg:   Tolerance for pcg, positive scalar < 1, 
            %               default 1e-3 and adpatively update.
            %   @tol_pcg_up:Factor to update tol_pcg, positive scalar < 1, 
            %               default 0.7.
            %   @maxiter:   Maximal iteration, positive scalar, default 10000.
            %   @disp_freq: Frequency for information display, default 1, 
            %               set to "inf" if you don't want any display.
            %   @showfig:   Showing or not showing figures of difference of 
            %               Theta per disp_freq, default false.

            arguments (Input)
                opts.Y double
                opts.lamb1 (1, 1) {mustBeNumeric, mustBeNonnegative} = 0.1
                opts.lamb2 (1, 1) {mustBeNumeric, mustBeNonnegative} = 100
                opts.lamb3 (1, 1) {mustBeNumeric, mustBePositive} = 0.5
                opts.beta_ (1, 1) {mustBeNumeric, mustBePositive} = 1
                opts.epsilon (1, 1) {mustBeNumeric, mustBeNonnegative} = 0.01
                opts.tol (1, 1) {mustBeNumeric, mustBePositive} = 1e-2
                opts.tol_pcg (1, 1) {mustBeNumeric, mustBePositive} = 1e-3
                opts.tol_pcg_up (1, 1) {mustBeNumeric, mustBePositive} = 0.7
                opts.maxiter (1, 1) {mustBeNumeric, mustBePositive} = 5000
                opts.disp_freq (1, 1) {mustBeNumeric} = 1
                opts.showfig (1, 1) = false
            end

            % Data.
            obj.Y = opts.Y;
            % Parameters and options.
            obj.lamb1 = opts.lamb1;
            obj.lamb2 = opts.lamb2;
            obj.lamb3 = opts.lamb3;
            obj.epsilon = opts.epsilon;
            obj.beta_ = opts.beta_;
            obj.tol = opts.tol;
            obj.tol_pcg = opts.tol_pcg;
            obj.tol_pcg_up = opts.tol_pcg_up;
            obj.maxiter = opts.maxiter;
            obj.disp_freq = opts.disp_freq;
            obj.showfig = opts.showfig;

            % Obtain dimension and sample size.
            [obj.T, obj.d] = size(obj.Y);
            obj.Ym = reshape(obj.Y', 1, obj.d, obj.T);
            obj.Ymt = permute(obj.Ym, [2, 1, 3]);

            % Intermediate variables for acceleration.
            obj.first_col(obj.d * obj.d * obj.T, 1) = 0;
            obj.third_col(obj.d * obj.d * obj.T, 1) = 0;

            % Helper variable indicating indices of diagonal elements.
            obj.idx = (1 : obj.d+1 : obj.d^2).' + obj.d^2 .* (0:obj.T-1);

            % Compute Sigma.
            Y_reshaped = reshape(obj.Y, [1, obj.T, obj.d]);
            obj.Sigma = bsxfun(@times, permute(Y_reshaped, [3, 1, 2]), ...
                permute(Y_reshaped, [1, 3, 2]));
            obj.Sigma = (obj.Sigma + permute(obj.Sigma, [2, 1, 3])) / 2;

            % Compute SqrtSigma.
            obj.SqrtSigma = GFDtL.SqrtMat(obj.Sigma);

            % Initialization.
            obj.Initialization();

            % Compute initilized objective value.
            obj.PrimObjVal();
            obj.GetZeta_();
            obj.DualObjVal();
        end

        % The function to call the algorithm to solve the problem.
        run(obj)

        % The function to compute the regularizer in the objective function in
        % the primal problem.
        function val = R(obj, x)
            val = x .* (x <= obj.lamb3) + ...
                (x.^2 - obj.lamb3^2 + obj.lamb3) .* (x > obj.lamb3);
        end
        
        % The function to compute \mathcal{G} that appears in the objective
        % function in the dual problem.
        function val = G(obj, x)
            val = min(-max(x - obj.lamb2, 0) .* obj.lamb3, obj.lamb2 .* ...
                ((max(obj.lamb3 - x/2/obj.lamb2, 0)).^2 - x.^2/(2*obj.lamb2).^2 ...
                - obj.lamb3^2 + obj.lamb3));
        end
        
        % The function to initialize the primal / dual variables.
        Initialization(obj)
        
        % The function to compute the tridiagonal system in Theta update.
        RHS = TriSys(obj, x)

        % The function to compute the loss function.
        val = lossfunc(obj, Y)

        % The function to compute the primal objective value.
        PrimObjVal(obj)

        % The function to compute the dual objective value.
        DualObjVal(obj)

        % The function to update Theta.
        Theta_Up(obj)

        % The function to update V.
        V_Up(obj)

        % The function to update Upsilon.
        Upsilon_Up(obj)

        % The function to update D.
        D_Up(obj)
        
        % The function to obtain auxiliary variable zeta_.
        GetZeta_(obj)

        % The function to compute dual infeasibility.
        dfeas = DualFeas(obj)

        % Plot.
        plot(obj, figID, time_)

        function GetICs(obj)
            %% Compute AIC / BIC for the estimator.
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
            %   Est.GetBIC;
            %   Est.AIC
            %   Est.BIC
            %
            % Normally, this function needn't be called by the user as it
            % will be called after the function `run' is finished and the
            % problem is feasible.
            % If the function `run' is interupted by the user, then
            % `GetICs' should be called to compute the AIC / BIC.

            diff_Theta = diff(obj.Theta_k, 1, 3);
            threshold_ = 1e-6;

            idx_diff = ~repmat(eye(obj.d), [1, 1, obj.T - 1]);
            idx_1 = 1:obj.d^2;
            K = sum(abs(diff_Theta(idx_diff)) >= threshold_) ...
                + sum(abs(obj.Theta_k(idx_1(~eye(obj.d)))) >= threshold_, 'all');
            obj.AIC = 2 * K + 2 * obj.lossval;
            obj.BIC = K * log(obj.T) + obj.lossval;
        end

        function breaks = EstBreaks(obj)
            %% Use adpative threshold to estimate the breaks.
            % The breaks are obtained by picking those time slots whose D_k
            % are nonzero.
            %
            % - Usage:
            %   breaks = Est.EstBreaks;
            %   breaks % OR 
            %   Est.breaks

            if obj.infeas
                breaks = (1:obj.T)';
            else
                tmp_norm = vecnorm(reshape(obj.D_k, [], obj.T-1));
                breaks = find(tmp_norm >= 1e-6)';
                breaks = mergeBreaks(breaks);
            end

            obj.breaks = breaks;
        end
    end

    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% The followings are declarations of helper functions as      %%%%
        %%%% static functions.                                           %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Xvec = vec(X)
            %% Given a high dimensional array X, return the vectorized column vector.
            % 
            % - Input:
            %   @X:     A high dimensional array.
            %
            % - Output:
            %   @Xvec:  A vectorized column vector.

            Xvec = X(:);
        end

        function Xsvec = svec(X)
            %% Given X containing several symmetric matrices, return the symmetric vectorized vectors.
            %
            % - Input:
            %   @X:     A high dimensional array of size [d, d, T], for
            %           each t, X(:, :, t) is a symmetric matrix.
            %
            % - Output:
            %   @Xsvec: A matrix of size [d(d+1)/2, T].

            [d, ~, T] = size(X);
            idx_multi_tril = repmat(~triu(ones(d, d), 1), [1, 1, T]);
            Xsvec = reshape(X(idx_multi_tril), [], T);
        end
        
        function SqrtM = SqrtMat(M)
            %% Given a positive semidefinite matrix M, compute M^{1/2}.
            % Compare to using the built-in function "sqrtm", we use "eig" to 
            % compute this matrix square root to avoid unnecessary mistakes.
            % This function also accept M as a tensor, i.e., M = [M(i)] where 
            % M(i) is a matrix for each i.
            % To achieve this compatibility, we shall use "pageeig".
            %
            % - Input:
            %   @M:     A matrix or a three-dimension tensor.
            %
            % - Output:
            %   @SqrtM: A matrix M^{1/2} or a three-dimension tensor consists of
            %           matrix square root.
            
            % Compute eigendecomposition.
            [Vs, Ds] = pageeig(M);
            % Drop some numerical error terms.
            Ds(abs(Ds) <= 1e-10) = 0;
            % Compute matrix sqrt.
            SqrtM = pagemtimes(pagemtimes(Vs, sqrt(Ds)), permute(Vs, [2, 1, 3]));
        end

        function Xinv = Inv(X)
            %% Compute the inverse of a matrix X via eigendecomposition.
            % Compare to the built-in inv, this function first symmetrize the 
            % matrix. Then an eigendecomposition is applied to compute the 
            % inverse.
            % - Input:
            %   @X:     The input matrix, suppose to be symmetric or near 
            %           symmetric.
            % - Output:
            %   @Xinv:  The inverse of X.

            X = (X + X') / 2;
            [V, D] = eig(X);
            Xinv = V * diag(1 ./ diag(D)) * V';
            Xinv = (Xinv + Xinv') / 2;
        end
    end
end