function [evalstot, modestot, dynamicstot, recons] = ...
    model_builder(gin, dmd_thresh, dt, num_pred)

    global num_cores

    [num_rows, num_cols, num_comps] = size(gin);
    
    evalstot = {};
    modestot = {};
    dynamicstot = {};
    recons = zeros(num_rows, num_pred, num_comps);

    fprintf("   Number of components: %d\n", num_comps)
    for ll=1:num_comps
        % Center data for this connected component
        gin_mean = mean(mean(gin(:, :, ll)));
        gin(:, :, ll) = gin(:, :, ll) - gin_mean;

        % DMD step
        gm = gin(:, 1:end-1, ll);
        gp = gin(:, 2:end, ll);
        [u, s, v] = svd(gm, 'econ');
        sr = diag(s);
        indskp = log10(sr / max(sr)) > dmd_thresh;
        srd = sr(indskp);
        kmat = gp * v(:, indskp) * diag(1./srd) * u(:, indskp)';
        [evecs, evals] = eig(kmat);
        phi = evecs\gm;

        % Reconstruct training data (1-step recon)
        recons(:, 1, ll) = real(evecs * phi(:, 1));
        recons(:, 2:num_cols, ll) = real(evecs * evals * phi);
        
        % Forecast from last profile to end
        gint = phi(:, end);
        num_steps = num_pred-num_cols;
        if num_cores>0
            evals_list = zeros(num_rows, num_steps);
            for ii=1:num_steps
                evals_list(:, ii) = diag(evals).^ii;
            end
            parfor jj=1:num_steps
                recons(:, num_cols+jj, ll) = real(evecs * diag(evals_list(:, jj)) * gint);
            end
        else
            evls_cur = evals.^2;
            for jj = num_cols:num_pred
                recons(:, jj, ll) = real(evecs * evls_cur * gint);
                evls_cur = evls_cur.*evals;
            end
        end

        % Reintroduce the means
        recons = recons + gin_mean;
        gin(:, :, ll) = gin(:, :, ll) - gin_mean;

        % Save all eigenvales/modes/functions
        evalstot{end+1} = diag(evals);
        modestot{end+1} = evecs;
        dynamicstot{end+1} = phi;

        fprintf("   Component %d\n", ll)
    end
    
end
