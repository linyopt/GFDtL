function DualObjVal(obj)
    %% Compute the dual objective value.

    % Three parts:
    % 1. -\tr(Wt Wt^T) / 2 with Wt = (Xt Xt^T)^{1/2} Theta_t^k.
    % 2. \epsilon \tr(zeta_). Here we compute trace by summing its diagonal
    %    elements.
    % 3. G part.
    obj.dualval = -sum(pagemtimes(obj.Theta_k, obj.Theta_k) .* obj.Sigma, 'all') / 2 ...
        + obj.epsilon * sum(obj.zeta_(obj.idx), 'all') ...
        + sum(obj.G(vecnorm(reshape(obj.Z_k, [], obj.T - 1)))); 
end