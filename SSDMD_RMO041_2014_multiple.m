%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Scale separated DMD on the RO041 2014 dataset
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

% Set up for parallel processing
if isempty(gcp('nocreate'))
    parpool(8); % Be sure to set number of cpu cores here
end

% Load sounder data
station = 'RO041';
file_date = '2014-01-01(365)';
data = load(strcat(station, '_Digisonde_', file_date));
h_low = 1;                     
h_high = 350;
ne_full = data.ne(h_low:h_high, :);
heights = data.heights(h_low:h_high);
full_times = data.times;

% Load IRI
iri = load(strcat(station, '_IRI_', file_date));
iri_hmin = find(iri.alt_km==heights(1));
iri_hmax = find(iri.alt_km==heights(end));
ne_iri = iri.ne(iri_hmin:iri_hmax, :);
fof2_iri = iri.foF2;
hmf2_iri = iri.hmF2;

% Load Didbase profile characteristics
chars = load(strcat(station, '_FastChars_', file_date));
fof2_data = chars.fof2;
hmf2_data = chars.hmf2;

% Some parameters for the SSDMD model
day = 96;
dt = 1./day;
wave_type = 'coif4';
dmd_tol = -6.0;
corr_tol = -1.75;
nd_train = 10;
nd_test = 2;
nd = nd_train + nd_test;
wave_levels = floor(log2(nd_train*day + 1));
% Random initial conditions to train on during the year
bad1 = [28978, 30358];
bad2 = [16911, 19307];
n_rstart = 30;
rng(1997, 'twister');
s = RandStream('mlfg6331_64');
r1 = datasample(s, 3000:17300, 15, 'Replace', false);
r2 = datasample(s, 24700:29600, 5, 'Replace', false);
r3 = datasample(s, 30350:33800, 10, 'Replace', false);
rstart = cat(2, r1, r2, r3);


%% Build/test SSDMD models for random initial starts over the year
n_comps = zeros(n_rstart, 1);
quick_mae = zeros(n_rstart, 1);
full_fof2_ssdmd = zeros(n_rstart, day*nd);
full_fof2_data = zeros(n_rstart, day*nd);
full_fof2_iri = zeros(n_rstart, day*nd);
full_hmf2_ssdmd = zeros(n_rstart, day*nd);
full_hmf2_data = zeros(n_rstart, day*nd);
full_hmf2_iri = zeros(n_rstart, day*nd);
full_test_times = zeros(n_rstart, day*nd);
t0=tic;
parfor ii=1:n_rstart
    % Set indices in the larger dataset for training and test subsets
    istart = rstart(ii);
    start_ix = istart;
    stop_ix = start_ix + nd*day - 1;
    fprintf("Start date: %s\n", full_times(istart,:))
    
    ne = ne_full(:, start_ix:stop_ix);
    [num_rows, num_cols] = size(ne);
    
    % Split sounder train/test data
    train_stop_ix = nd_train*day;
    train_data = ne(:, 1:train_stop_ix);
    test_start_ix = train_stop_ix + 1;
    test_stop_ix = train_stop_ix + nd_test*day;
    
    % Fit SSDMD model and make forecast over test data range
    model = parssdmd(train_data, wave_levels, wave_type, dmd_tol, ...
        corr_tol, day, dt, num_cols, heights);
    
    % Save forecasted parameters for this period
    full_fof2_ssdmd(ii, :) = model.fof2;
    full_hmf2_ssdmd(ii, :) = model.hmf2;
    full_test_times(ii, :) = iri.times(start_ix:stop_ix);

    % Save the test data for this period
    full_fof2_data(ii, :) = fof2_data(start_ix:stop_ix);
    full_hmf2_data(ii, :) = hmf2_data(start_ix:stop_ix);

    % Save the IRI data for this period
    full_fof2_iri(ii, :) = fof2_iri(start_ix:stop_ix);
    full_hmf2_iri(ii, :) = hmf2_iri(start_ix:stop_ix);

    n_comps(ii) = model.num_comps;
    quick_mae(ii) = mean(abs(full_hmf2_data(ii,:) - full_hmf2_ssdmd(ii,:)), 'omitnan')
end
fprintf('[COMPLETE] '); toc(t0)
fprintf('\n avg. num_comps: %d\n', mean(n_comps))
fprintf('\n avg. mae: %d\n', mean(quick_mae))


%%
save('Run_RO041_2014')

