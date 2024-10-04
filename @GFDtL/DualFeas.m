function dfeas = DualFeas(obj)
    %% Compute dual infeasibility.
    
    % Compute eigenvalues of zeta_.
    eigs_zeta_ = reshape(pageeig(obj.zeta_), obj.d, obj.T);
    % Find the absolute value of the minimum negative eigenvalues.
    min_eigs = abs(min(min(eigs_zeta_), 0));
    % Compute the norm of zeta_.
    norms_zeta_ = vecnorm(reshape(obj.zeta_, [], obj.T));
    % Compute dfeas1.
    dfeas1 = max(min_eigs ./ (norms_zeta_ + 1));

    % Compute maximum value of Y_k.
    max_Y_k = max(abs(obj.Y_k), [], 'all');
    % Compute dfeas2.
    dfeas2 = max(max_Y_k - obj.lamb1, 0) / (1 + max_Y_k);

    % Compute dfeas.
    dfeas = max([dfeas1, dfeas2]);
end