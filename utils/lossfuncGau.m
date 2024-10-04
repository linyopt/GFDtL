function lossfunc = lossfuncGau(Theta, Y)
%% Compute the loss function value of GFGL.
% The input Y can be a series of data, then the output will be the mean of
% the loss function values.
%
% - Usage:
%   lossfunc = lossfuncGau(Theta, Y);
%
% - Input:
%   @Theta:     The estimated Theta.
%   @Y:         The data.
%
% - Output:
%   @lossfunc:  The loss function value.

% Dimension and sample size.
[d, ~, T] = size(Theta);

% The number of datasets in Y.
num_Y = size(Y, 3);

% The loss function values.
vals(num_Y) = 0;

% For each dataset in Y.
for tt = 1:num_Y
    % Compute the Sigma.
    Y_reshaped = reshape(Y(:, :, tt), [1, T, d]);
    S = bsxfun(@times, permute(Y_reshaped, [3, 1, 2]), permute(Y_reshaped, [1, 3, 2]));

    % Compute the corresponding loss function value.
    % The following four lines come from the codes of the GFGL paper.
    lossfunc=-log(det(squeeze(Theta(:,:,1))))+trace(squeeze(S(:,:,1))*squeeze(Theta(:,:,1)));
    for t=2:T
        lossfunc=lossfunc-log(det(squeeze(Theta(:,:,t))))+trace(squeeze(S(:,:,t))*squeeze(Theta(:,:,t)));
    end

    vals(tt) = lossfunc;
end
% Compute mean and return.
lossfunc = mean(vals);
end

