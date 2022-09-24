%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Scale separated DMD on the GA762 2022 dataset
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

% Load sounder data
station = 'GA762';
file_date = '2022-06-01(015)';
data = load(strcat(station, '_Digisonde_', file_date));
h_low = 1;                     
h_high = 350;
ne_full = data.ne(h_low:h_high, :);
heights = data.heights(h_low:h_high);
full_times = data.times;

% Load IRI
iri = load(strcat(station, '_IRI_', file_date));
ne_iri = iri.ne(151:500, :);
fof2_iri = iri.foF2;
hmf2_iri = iri.hmF2;

% Load Didbase profile characteristics
chars = load(strcat(station, '_FastChars_', file_date));
fof2_data = chars.fof2;
hmf2_data = chars.hmf2;

% Some parameters for the SSDMD model
day = 192;
dt = 1./day;
wave_type = 'coif4';
dmd_tol = -6.0;
corr_tol = -1.95;
nd_train = 10;
nd_test = 2;
nd = nd_train + nd_test;
wave_levels = wmaxlev(nd*day, wave_type);
start_ix = 1;
stop_ix = nd*day;

% Run model
ne = ne_full(:, start_ix:stop_ix);
[num_rows, num_cols] = size(ne);

% Split sounder train/test data
train_stop_ix = nd_train*day;
train_data = ne(:, 1:train_stop_ix);
test_start_ix = train_stop_ix + 1;
test_stop_ix = train_stop_ix + nd_test*day;

% Fit SSDMD model and make forecast over test data range
model = ssdmd(train_data, wave_levels, wave_type, dmd_tol, ...
    corr_tol, day, dt, num_cols, heights);

% Save forecasted parameters for this period
fof2_ssdmd = model.fof2(test_start_ix:test_stop_ix);
hmf2_ssdmd = model.hmf2(test_start_ix:test_stop_ix);

% Set the test data for this period
fof2_data = fof2_data(test_start_ix:test_stop_ix);
hmf2_data = hmf2_data(test_start_ix:test_stop_ix);

% Set the IRI data for this period
fof2_iri = fof2_iri(test_start_ix:test_stop_ix);
hmf2_iri = hmf2_iri(test_start_ix:test_stop_ix);

% Compute MAEs
fof2_ssdmd_mae = mean(abs(fof2_data - fof2_ssdmd), 'omitnan');
hmf2_ssdmd_mae = mean(abs(hmf2_data - hmf2_ssdmd), 'omitnan');
fof2_iri_mae = mean(abs(fof2_data - fof2_iri), 'omitnan');
hmf2_iri_mae = mean(abs(hmf2_data - hmf2_iri), 'omitnan');



%% Plot forecast for station
ms = 20;
lw = 6;
fs = 24;

% fof2
figure
tl = tiledlayout(2, 1, 'TileSpacing','tight', 'Padding', 'tight');
nexttile;
hold on
plot(fof2_data, 'r.', 'MarkerSize', ms)
plot(fof2_ssdmd, 'b-', 'LineWidth', lw)
plot(fof2_iri, 'k-', 'LineWidth', lw)
grid on;
hold off
labelssdmd = sprintf('SSDMD (MAE: %.2f)', fof2_ssdmd_mae);
labeliri = sprintf('IRI (MAE: %.2f)', fof2_iri_mae);
legend('Measured', labelssdmd, labeliri, 'Location','northwest')
ylabel('foF2 (MHz)', 'FontWeight', 'bold')
xticks([0.5*day:day:nd_test*day-0.5*day])
xticklabels([datestr(makeDate(iri.times(test_start_ix:day:end)), 'yyyy/mm/dd')])
h = gca;
h.FontSize = fs;
% hmf2
nexttile;
hold on
plot(hmf2_data, 'r.', 'MarkerSize', ms)
plot(hmf2_ssdmd, 'b-', 'LineWidth', lw)
plot(hmf2_iri, 'k-', 'LineWidth', lw)
grid on;
hold off
labelssdmd = sprintf('SSDMD (MAE: %.2f)', hmf2_ssdmd_mae);
labeliri = sprintf('IRI (MAE: %.2f)', hmf2_iri_mae);
legend('Measured', labelssdmd, labeliri, 'Location','northwest')
ylabel('hmF2 (km)', 'FontWeight', 'bold')
xticks([0.5*day:day:nd_test*day-0.5*day])
xticklabels([datestr(makeDate(iri.times(test_start_ix:day:end)), 'yyyy/mm/dd')])
h = gca;
h.FontSize = fs;


%% Total reconstructions with hmf2 param
ms = 20;
lw = 6;
fs = 24;

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
xl = xline(nd_train*day, 'w--', 'linewidth', 4);
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
xl = xline(nd_train*day, 'w--', 'linewidth', 4);
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



