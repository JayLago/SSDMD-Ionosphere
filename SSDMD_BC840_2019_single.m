%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Scale separated DMD on the BC840 2019 dataset for a short snippet
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

% Load sounder data
station = 'BC840';
file_date = '2019-01-01(365)';
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
day = 288;
dt = 1./day;
wave_type = 'coif4';
dmd_tol = -6.0;
corr_tol = -1.95;
nd_train = 10;
nd_test = 2;
nd = nd_train + nd_test;
num_days = nd;
wave_levels = floor(log2(nd*day));
start_ix = 79800;   % October 5 2019
stop_ix = start_ix+nd*day;

% Get test data from large data set
ne = ne_full(:, start_ix:stop_ix);
[num_rows, num_cols] = size(ne);
ne_iri = ne_iri(:, start_ix:stop_ix);
fof2_data = fof2_data(start_ix:stop_ix);
hmf2_data = hmf2_data(start_ix:stop_ix);

% Split sounder train/test data
train_stop_ix = nd_train*day;
test_stop_ix = nd*day;
train_data = ne(:, 1:train_stop_ix);

% Fit SSDMD model and make forecast over test data range
model = ssdmd(train_data, wave_levels, wave_type, dmd_tol, ...
    corr_tol, day, dt, num_cols, heights);



%% Total reconstructions with hmf2 param plotted on top
fs = 36;

% Shift hmf2 so the indices line up with data in plot
hmf2_ = hmf2_data - double(heights(1));

tiledlayout(2, 1, 'TileSpacing','Compact', 'Padding', 'Compact');
nexttile;
s1 = surf(ne);
hold on;
shading flat;
view([0, 90]);
z_max = max(max(ne)+1)*ones(num_cols);
x_line = linspace(1, num_cols, num_cols);
hl = line(x_line, hmf2_, z_max, 'Color', 'k', 'LineWidth', 2);
c = colorbar;
colormap('jet')
c.Label.String = 'Plasma Frequency (MHz)';
caxis([0 6]);
xl = xline(nd_train*day, 'm--', 'linewidth', 4);
title('Measurement', 'FontSize', fs);
h = set(gca,'FontSize', fs);
set(h,'Interpreter','LaTeX')
xticks([1:day:num_cols])
xticklabels([])
xlim([0, num_cols])
ylim([0, num_rows])
xlabel([])
yticks([1:100:num_rows])
yticklabels([heights(1:100:end)])
ylabel('Height (km)', 'FontSize', fs)
legend([hl], 'hmF2')
hold off


% Shift hmf2 so the indices line up with data in plot
hmf2_ = model.hmf2 - double(heights(1));

nexttile;
s2 = surf(model.full_recon);
hold on;
shading flat;
view([0, 90]);
z_max = max(max(ne)+1)*ones(num_cols);
x_line = linspace(1, num_cols, num_cols);
hl = line(x_line, hmf2_, z_max, 'Color', 'k', 'LineWidth', 2);
c = colorbar;
colormap('jet')
c.Label.String = 'Plasma Frequency (MHz)';
caxis([0 6]);
xl = xline(nd_train*day, 'm--', 'linewidth', 4);
title('SSDMD', 'FontSize', fs);
h = set(gca,'FontSize', fs);
set(h,'Interpreter','LaTeX')
xticks([1:day:num_cols])
xticklabels([0:num_days])
xlim([0, num_cols])
ylim([0, num_rows])
xlabel('Time (days)', 'FontSize', fs)
yticks([1:100:num_rows])
yticklabels([heights(1:100:end)])
ylabel('Height (km)', 'FontSize', fs)
legend([hl], 'hmF2')
hold off;


