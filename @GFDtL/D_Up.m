function D_Up(obj)
    %% Update D_k.
    
    % Compute Xi.
    Xi = diff(obj.Theta_k, 1, 3) - obj.Z_k ./ obj.beta_;
    % Compute norm of Xi.
    Norms_Xi = vecnorm(reshape(Xi, [], obj.T - 1));
    
    % Compute iota1 and iota2; see the paper for details.
    iota1 = min(max(Norms_Xi - obj.lamb2 / obj.beta_, 0), obj.lamb3);
    iota2 = max(Norms_Xi ./ (1 + 2 * obj.lamb2 / obj.beta_), obj.lamb3);
    % Compute the two corresponding values; see the paper for details.
    Vals1 = 1/2.*(iota1 - Norms_Xi).^2 + obj.lamb2 / obj.beta_ .* iota1;
    Vals2 = 1/2 .* (iota2 - Norms_Xi).^2 + obj.lamb2 / obj.beta_ .* (iota2.^2 - obj.lamb3^2 + obj.lamb3);
    
    % Declare the updated D_k.
    D_kp1 = zeros(obj.d, obj.d, obj.T - 1);
    % Find the indices of zero Xi's.
    % This is to handle some edge cases.
    idx0 = Norms_Xi == 0;
    % Find the indices for the first cases.
    idx = Vals1 <= Vals2;

    % Compute iota / norm, i.e., the coefficients of Xi.
    % Here we only consider nonzero Xi's and use idx to select the
    % corresponding parts.
    % For example, idx refers to the first case, while ~idx refers to the
    % second case.
    temp1 = iota1(idx & ~idx0) ./ Norms_Xi(idx & ~idx0);
    temp2 = iota1(~idx & ~idx0) ./ Norms_Xi(~idx & ~idx0);
    
    % Compute the real D_kp1 by multiplying the coefficients and Xi.
    D_kp1(:, :, idx & ~idx0) = bsxfun(@times, Xi(:, :, idx & ~idx0), reshape(temp1, 1, 1, []));
    D_kp1(:, :, ~idx & ~idx0) = bsxfun(@times, Xi(:, :, ~idx & ~idx0), reshape(temp2, 1, 1, []));
    
    % Symmetrization.
    obj.D_k = (D_kp1 + permute(D_kp1, [2, 1, 3])) / 2;
end