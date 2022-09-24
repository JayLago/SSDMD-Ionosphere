function [model] = ssdmd(train_data, wave_levels, wave_type, ...
    dmd_tol, corr_tol, cols_per_day, dt, num_pred, heights)

    % Center the data
    meanNe = mean(train_data, 2);
    train_data = train_data - meanNe;
    
    % Separate scales and compute correlation values in training data
    fprintf('Separating scales...\n')
    t0 = tic;
    [sep_scales, corr_vals] = separateScales(train_data, wave_levels, wave_type);
    fprintf('[COMPLETE] '); toc(t0)
    
    
    % Group the correlated scales together to get connected components
    fprintf('Correlating scales...\n')
    t0 = tic;
    [comps, num_comps, graph] = getConnectedComps(corr_tol, sep_scales, corr_vals, train_data);
    fprintf('[COMPLETE] '); toc(t0)
    
    
    % Find the average signal and standard deviations of the noise from each wavelet group
    fprintf('Computing averages and flucuations...\n')
    t0 = tic;
    [gavgs, gstds] = separateAverages(comps, cols_per_day);
    fprintf('[COMPLETE] '); toc(t0)
    
    
    % Build DMD model for averages over 24 hour windows of time
    fprintf('Building DMD model for averages...\n')
    t0 = tic;
    [evals, modes, dynamics, recons] = runDMD(gavgs, dmd_tol, dt, num_pred);
    fprintf('[COMPLETE] '); toc(t0)
    
    % Reconstruction
    full_recon = sum(recons, 3) + meanNe;
    
    % Compute foF2 and hmF2 parameters
    [fof2, hmf2_idx] = max(full_recon, [], 1);
    hmf2 = cast(heights(hmf2_idx), 'double');
    
    % Populate output struct
    model.sep_scales = sep_scales;
    model.corr_vals = corr_vals;
    model.num_comps = num_comps;
    model.comps = comps;
    model.graph = graph;
    model.avgs = gavgs;
    model.flucs = gstds;
    model.evals = evals;
    model.modes = modes;
    model.dynamics = dynamics;
    model.recons = recons;
    model.full_recon = full_recon;
    model.meanNe = meanNe;
    model.fof2 = fof2;
    model.hmf2 = hmf2;

end
















