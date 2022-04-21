function [model] = ssdmd(train_data, wave_levels, wave_type, ...
    dmd_tol, corr_tol, cols_per_day, dt, num_pred)

[num_rows, num_cols] = size(train_data);


%% Separate scales and compute correlation values in training data
fprintf('Separating scales...\n')
tic;
[sep_scales, corr_vals] = correlation_mat_maker(train_data, wave_levels, wave_type);
fprintf('[COMPLETE] ')
toc


%% Generate wavelet separation of the total signal 
fprintf('Correlating scales...\n')
tic
[comps, num_comps, graph] = wavelet_separation(corr_tol, sep_scales, corr_vals, train_data);
fprintf('[COMPLETE] ')
toc


%% Find the average signal and standard deviations of the noise from each wavelet group
fprintf('Computing averages and flucuations...\n')
gavgs = zeros(size(comps));
gflucs = zeros(size(comps));
gstd = zeros(num_rows, cols_per_day, num_comps);
tdays = num_cols/cols_per_day;
tic
for ll = 1:num_comps
    for jj = 1:tdays
        gavgs(:, :, ll) = gavgs(:,:,ll) + ...
            repmat(comps(:, 1+(jj-1)*cols_per_day:jj*cols_per_day,ll), 1, tdays) / tdays;
    end
    gflucs(:, :, ll) = comps(:, :, ll) - gavgs(:, :, ll);
    gvar = zeros(num_rows, cols_per_day);
    for jj=1:tdays
        gvar = gvar + gflucs(:, 1+(jj-1)*cols_per_day:jj*cols_per_day, ll).^2;
    end    
    gstd(:, :, ll) = sqrt(gvar/(tdays-1));
end
fprintf('[COMPLETE] ')
toc


%% Build DMD model for averages over 24 hour windows of time
fprintf('Building DMD model for averages...\n')
tic;
[evals, modes, dynamics, recons] = model_builder(gavgs, dmd_tol, dt, num_pred);
fprintf('[COMPLETE] ')
toc


%% Populate output struct
model.sep_scales = sep_scales;
model.corr_vals = corr_vals;
model.num_comps = num_comps;
model.comps = comps;
model.graph = graph;
model.avgs = gavgs;
model.flucs = gstd;

model.evals = evals;
model.modes = modes;
model.dynamics = dynamics;
model.recons = recons;
model.full_recon = sum(model.recons, 3);

end
















