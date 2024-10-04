function Theta_Up(obj)
%% Theta update in nonprox version, solved inexactly by pcg.

% Compute Z_diff = Z_t - Z_{t-1}, where the first one and the last one are
% different.
Z_diff = zeros(obj.d, obj.d, obj.T);
% The first one.
Z_diff(:, :, 1) = obj.Z_k(:, :, 1);
% The middle ones.
Z_diff(:, :, 2:obj.T-1) = diff(obj.Z_k, 1, 3);
% The last one.
Z_diff(:, :, obj.T) = - obj.Z_k(:, :, obj.T-1);

% Compute D_diff = D_t - D_{t-1}.
D_diff = zeros(obj.d, obj.d, obj.T);
% The first one.
D_diff(:, :, 1) = obj.D_k(:, :, 1);
% The middle ones.
D_diff(:, :, 2:obj.T-1) = diff(obj.D_k, 1, 3);
% The last one.
D_diff(:, :, obj.T) = - obj.D_k(:, :, obj.T-1);

% Compute Gamma.
Gamma_k = repmat(eye(obj.d), [1, 1, obj.T]) + obj.A_k + obj.Y_k - Z_diff ...
    + obj.beta_ * obj.V_k + obj.beta_ * obj.Upsilon_k - obj.beta_ * D_diff;

% Call pcg.
[Theta_k_vec, ~] = pcg(@(x)obj.TriSys(x), obj.vec(Gamma_k), obj.tol_pcg, obj.d^2, [], [], obj.vec(obj.Theta_k));

% Reshape into a sequence of matrices.
Theta_kp1 = reshape(Theta_k_vec, obj.d, obj.d, obj.T);

% Symmetrization.
obj.Theta_k = (Theta_kp1 + permute(Theta_kp1, [2, 1, 3])) / 2;
end