function PrimObjVal(obj)
    %% Compute the primal objective value.

    % Compute norm of first-order difference of Theta.
    obj.norm_diff = vecnorm(reshape(diff(obj.Theta_k, 1, 3), [], obj.T-1));

    % Check infeasibility.
    if max(obj.norm_diff) >= obj.lamb3
        obj.infeas = true;
    end

    % Compute lossval.
    % Here we compute the trace by summing the diagonal elements.
    obj.lossval = sum(pagemtimes(obj.Theta_k, obj.Theta_k) .* obj.Sigma, 'all') / 2 - ...
            sum(obj.Theta_k(obj.idx), 'all');
    
    % Compute primal value.
    % 1. lossval.
    % 2. Off-diagonal l1.
    % 3. Gourp fused lasso.
    obj.primval = obj.lossval + obj.lamb1 * (sum(abs(obj.Theta_k), 'all') - sum(abs(obj.Theta_k(obj.idx)), 'all')) + ...
            obj.lamb2 * sum(obj.R(obj.norm_diff));
end
