function Initialization(obj)
    %% Initilization of the algorithm.
    % This initialization method initializes the algorithm with inverse of 
    % sum of Sigma.
    
    % Compute sum of Sigma for the initialization of the algorithm.
    SumSigma = sum(obj.Sigma, 3);
    % Compute the inverse of sum of Sigma, it exists due to the
    % assumption that SumSigma is positive definite.
    InvSumSigma = obj.Inv(SumSigma);
    obj.Theta_k = repmat(InvSumSigma, [1, 1, obj.T]);
    obj.V_k = obj.Theta_k;
    obj.Upsilon_k = obj.Theta_k;
    obj.D_k = zeros(obj.d, obj.d, obj.T-1);
    obj.A_k = ones(obj.d, obj.d, obj.T);
    obj.Y_k = zeros(obj.d, obj.d, obj.T);
    obj.Z_k = zeros(obj.d, obj.d, obj.T - 1);
end