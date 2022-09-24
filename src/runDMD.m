function [evalstot, modestot, dynamicstot, recons] = runDMD(gin, dmd_thresh, dt, num_pred)

    [num_rows, num_cols, num_comps] = size(gin);
    
    evalstot = {};
    modestot = {};
    dynamicstot = {};
    recons = zeros(num_rows, num_pred, num_comps);

%     fprintf("\tNumber of components: %d\n", num_comps)
    for ll=1:num_comps
%         fprintf("\tComponent: %d\n", ll)
        gm = gin(:, 1:end-1, ll);
        gp = gin(:, 2:end, ll);
        [u, s, v] = svd(gm, 'econ');
        sr = diag(s);
        indskp = log10(sr / max(sr)) > dmd_thresh;
        srd = sr(indskp);
        kmat = gp * v(:, indskp) * diag(1./srd) * u(:, indskp)';
        [evecs, evals] = eig(kmat);
        phi = evecs\gm;
        recons(:, 1, ll) = real(evecs * phi(:, 1));
        recons(:, 2:num_cols, ll) = real(evecs * evals * phi);
        
        gint = phi(:, end);
        evls_cur = evals.^2;
        for jj = num_cols:num_pred
            recons(:, jj, ll) = real(evecs * evls_cur * gint);
            evls_cur = evls_cur.*evals;
        end

        evalstot{end+1} = diag(evals);
        modestot{end+1} = evecs;
        dynamicstot{end+1} = phi;
    end
    
end
