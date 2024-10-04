function RHS = TriSys(obj, x)
%% Compute the RHS value of a tridiagonal system that appears in Theta update.
% - Input:
%   @x:     The input for the tridiagonal system.
% - Output:
%   @RHS:   The RHS value of the tridiagonal system.

% As the input x is a vector (for pcg), we reshpae it into a matrix.
x_mat = reshape(x, obj.d, obj.d, obj.T);

% Compute the first part:
%   (XX^T \Theta + \Theta XX^T) / 2.
% In the code we use `Y` for data.
% Here we make use of Ymt and Ym, which are some reshaped Y and are defined
% in GFDtL.
% We use pagemtimes to compute XX^T \Theta.
prod_sum_Sig_xmat = pagemtimes(obj.Ymt, pagemtimes(obj.Ym, x_mat));
% Symmetrize it.
prod_sum_Sig_xmat = (prod_sum_Sig_xmat + permute(prod_sum_Sig_xmat, [2, 1, 3])) / 2;

% Compute the second part:
%   4\beta \Theta.
scaled_xmat = 4 * obj.beta_ * x_mat;
% The first one and the last one are 3\beta \Theta.
scaled_xmat(:, :, 1) = 3 * obj.beta_ * x_mat(:, :, 1);
scaled_xmat(:, :, obj.T) = 3 * obj.beta_ * x_mat(:, :, obj.T);
% All of them need to subtract the diagonal elements.
scaled_xmat(obj.idx) = scaled_xmat(obj.idx) - obj.beta_ * x_mat(obj.idx);
% Reshape into a vector.
second_col = reshape(prod_sum_Sig_xmat + scaled_xmat, [], 1);

% Both first_col and third_col are defined in GFDtL as a zero vector.
% Compute the first column:
%   -\beta \Theta_{t-1}.
obj.first_col(obj.d^2 + 1 : end) = -obj.beta_ * x(1 : end - obj.d^2);
% Compute the third column:
%   -\beta \Theta_{t+1}.
obj.third_col(1 : end - obj.d^2) = -obj.beta_ * x(obj.d^2 + 1 : end);

% Sum them to get RHS.
RHS = obj.first_col + second_col + obj.third_col;
end