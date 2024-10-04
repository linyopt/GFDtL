function GetZeta_(obj)
    %% Compute the auxiliary variable zeta_.

    % Declare zeta_.
    zeta_ = zeros(obj.d, obj.d, obj.T);
    % For the first one.
    zeta_(:, :, 1) = obj.Z_k(:, :, 1) - eye(obj.d) + obj.Sigma(:, :, 1) * ...
            obj.Theta_k(:, :, 1) - obj.Y_k(:, :, 1);
    % For the middle ones.
    zeta_(:, :, obj.T) = -obj.Z_k(:, :, obj.T-1) - eye(obj.d) ...
            + obj.Sigma(:, :, obj.T) * obj.Theta_k(:, :, obj.T) ...
            - obj.Y_k(:, :, obj.T);
    % For the last one.
    zeta_(:, :, 2:obj.T-1) = obj.Z_k(:, :, 2:obj.T-1) - obj.Z_k(:, :, 1:obj.T-2) - ...
        repmat(eye(obj.d), [1, 1, obj.T-2]) + pagemtimes(obj.Ymt(:, :, 2:obj.T-1), ...
        pagemtimes(obj.Ym(:, :, 2:obj.T-1), obj.Theta_k(:, :, 2:obj.T-1))) - ...
        obj.Y_k(:, :, 2:obj.T-1);
    
    % Symmetrization.
    obj.zeta_ = (zeta_ + permute(zeta_, [2, 1, 3])) / 2;
end