function run(obj)
    %% Run the main loop of algorithm.

    tic
    if ~isinf(obj.disp_freq)
        fprintf(1, "Iter  |    pobj   |   dobj    | rel_gap  |   dfeas  |p_rel_chg |d_rel_chg |\n");
    end

    iter = 1;
    while iter <= obj.maxiter
        % Theta update.
        Theta_km1 = obj.Theta_k;
        obj.Theta_Up();

        % V update.
        obj.V_Up();

        % Upsilon update.
        obj.Upsilon_Up();

        % D update.
        obj.D_Up();

        % Dual update.
        % A update.
        A_km1 = obj.A_k;
        obj.A_k = obj.A_k - 1.61 * obj.beta_ * (obj.Theta_k - obj.V_k);

        % Y update.
        Y_km1 = obj.Y_k;
        obj.Y_k = obj.Y_k - 1.61 * obj.beta_ * (obj.Theta_k - obj.Upsilon_k);
        obj.Y_k(obj.idx) = 0;

        % Z update.
        Z_km1 = obj.Z_k;
        obj.Z_k = obj.Z_k - 1.61 * obj.beta_ * (diff(obj.Theta_k, 1, 3) - obj.D_k);

        % Compute primal objective value.
        obj.PrimObjVal();

        % Compute dual infeasibility and dual objective value.
        % Comptue auxiliary variable zeta_.
        obj.GetZeta_();
        % Compute dual infeasibility.
        dfeas = obj.DualFeas();
        % Compute dual objective value.
        obj.DualObjVal();

        % Compute primal dual gap.
        gap = abs(obj.primval - obj.dualval) / ...
                (1 + abs(obj.primval) + abs(obj.dualval));

        % Compute relative successive changes.
        p_chg = norm(obj.vec(obj.Theta_k - Theta_km1)) / ...
            (1 + norm(obj.vec(obj.Theta_k)) + norm(obj.vec(Theta_km1)));
        A_chg = norm(obj.vec(obj.A_k - A_km1)) / ...
            (1 + norm(obj.vec(obj.A_k)) + norm(obj.vec(A_km1)));
        Y_chg = norm(obj.vec(obj.Y_k - Y_km1)) / ...
            (1 + norm(obj.vec(obj.Y_k)) + norm(obj.vec(Y_km1)));
        Z_chg = norm(obj.vec(obj.Z_k - Z_km1)) / ...
            (1 + norm(obj.vec(obj.Z_k)) + norm(obj.vec(Z_km1)));
        d_chg = max([A_chg, Y_chg, Z_chg]);

        % Information display.
        if mod(iter, obj.disp_freq) == 0
            fprintf(1, "%6.0f|%+5.4e|%+5.4e|%5.4e|%5.4e|%5.4e|%5.4e|\n", ...
                    iter, obj.primval, obj.dualval, gap, dfeas, p_chg, d_chg);
            if obj.showfig
                obj.plot
                drawnow
            end
        end

        % Check stop criterion.
        if (gap <= obj.tol && dfeas <= obj.tol) || ... normal stop
            (max([p_chg, d_chg]) <= obj.tol / 1000) % solution not changing

            if ~isinf(obj.disp_freq)
                fprintf("Terminate at iter %g. \n", iter);
            end
            break;
        end

        if mod(iter, 10) == 0 && obj.tol_pcg > 1e-12
            obj.tol_pcg = obj.tol_pcg * obj.tol_pcg_up;
        end

        if obj.infeas
            if ~isinf(obj.disp_freq)
                fprintf("Terminate at iter %g because of infeasibility. \n", iter);
            end
            return
        end

        iter = iter + 1;
    end
    if ~isinf(obj.disp_freq)
        fprintf(1, ' ======== Time: %s ======== \n\n\n', duration([0, 0, toc]));
    end

    obj.GetICs();
end