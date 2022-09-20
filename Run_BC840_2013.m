%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Scale separated DMD on the BC840 2013 dataset
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

% To run single-threaded set num cores to -1
global num_cores
num_cores = 8;
if isempty(gcp('nocreate'))
    parpool(num_cores);
end

saveFigs = false;
figDPI = 150;
fontSize = 32;
fontSizeBig = 42;

dataName = 'BC840_2013-10-15_2013-11-15';
load(dataName);
nh_clip = 250;
ne = ne(1:nh_clip, :);
heights = heights(1:nh_clip);

day = 96;           % 15 minute measurement cadence => 96 measurements per day
num_days = 25;      % Full dataset contains 25 days exactly
[num_rows, num_cols] = size(ne);
dt = 1./day;
num_pred = num_cols;

% Parameters of the scale-separated DMD
wave_type = 'coif4';
dmd_tol = -6.0;
corr_tol = -1.95;

% Split train/test data
num_days_train = 20;
train_start_day = 1;
train_stop_day = num_days_train*day;
gtrain = ne(:, train_start_day:train_stop_day);
gtest = ne(:, train_stop_day+1:end);
wave_levels = floor(log2(train_stop_day - train_start_day + 1));

% Run SSDMD on training data
model = ssdmd(gtrain, wave_levels, wave_type, dmd_tol, corr_tol, day, dt, num_pred);

%% Compute model errors
% SSDMD absolute percent error (pointwise)
abs_per_err_ssdmd = 100*abs((ne - model.full_recon)./ne);

% IRI absolute precent error  (pointwise)
iri = load('IRI_Ne_2013-10-15');
iri_ne = iri.ne;
abs_per_err_iri = 100*abs((ne - iri_ne)./ne);

% Simple median model absolute percent error (pointwise)
median_model = zeros(length(heights), day);
tmp_ne = ne(:, 1:train_stop_day);
for ii=1:day
    median_model(:, ii) = median(tmp_ne(:, ii:day:end), 2);
end
median_model = repmat(median_model, [1, num_days]);
abs_per_err_median = 100*abs((ne - median_model)./ne);


%% The full data set
tiledlayout(1, 1, 'Padding', 'Compact');
nexttile;
hold on
surf(ne);
shading flat
view([0, 90]);
c = colorbar;
c.Label.String = 'Plasma Frequency (MHz)';
c.FontSize = fontSize;
h = set(gca,'FontSize', fontSize);
xticks([1:day:num_cols])
xticklabels([0:num_days])
xlabel('Time (days)', 'FontSize', fontSize)
xlim([0, num_cols])
ylim([0, 250])
yticks([1:100:num_rows])
yticklabels([heights(1:100:end)])
ylabel('Height (km)', 'FontSize', fontSize)
title("Boulder, CO (2013/10/15 - 2013/11/09)", 'FontSize', fontSize);
hold off;

fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 10.7222 24.7778 7.8333])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_data_2013.eps', 'Resolution', figDPI)
    close(fig);
end


%% Hilbert spectrum
tmpne = ne - mean(mean(ne));
imf = emd(tmpne(230, :));
hht(imf, day, 'FrequencyLimits', [0 48])
h = set(gca,'FontSize', fontSize);
set(h,'Interpreter','LaTeX')
c = colorbar;
c.Label.String = 'Instantaneous Energy';
c.FontSize = fontSize;
% caxis([1, 10]);
xlabel("Time (days)", 'FontSize', fontSize)
ylabel("Frequency (cycles/day)", 'FontSize', fontSize)
title('Hilbert Spectrum');
ylim([-1 40])
fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 10.7639 16.8194 7.7917])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_data_hilbert_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end


%% All wavelet scale reconstructions
wave_scales = model.sep_scales;
[nrows, ncols, nscales] = size(wave_scales);
figure
tiledlayout(2, 3, 'TileSpacing','Compact', 'Padding', 'Compact');
hold on
for ii=1:6
    nexttile;
    surf(wave_scales(:, :, ii));
    shading flat
    view([0, 90]);
    h = set(gca,'FontSize', fontSizeBig);
    c = colorbar;
    yticks([1:100:num_rows])
    yticklabels([heights(1:100:end)])
    ylim([0, 250])
    ylabel('Height (km)', 'FontSize', fontSizeBig);
    xticks(1:3*day:ncols+1);
    xticklabels([0:3:num_days+1]);
    xlabel('Time (days)', 'FontSize', fontSizeBig);
    ylim([0, 250])
    xlim([0, train_stop_day])
    title(sprintf("Wavelet Scale %d", ii), 'FontSize', fontSizeBig);
end
hold off;
fig1 = gcf;
fig1.Units = 'inches';
set(fig1, 'position', [0 0 47.7778 18.5555])
if saveFigs
    exportgraphics(fig1, './results/Boulder/fig_wavelet_scale_reconstructions1_2013.eps', ...
        'Resolution', figDPI)
    close(fig1);
end

figure
tiledlayout(2, 3, 'TileSpacing','Compact', 'Padding', 'Compact');
hold on
for ii=7:nscales
    nexttile;
    surf(wave_scales(:, :, ii));
    shading flat
    view([0, 90]);
    h = set(gca,'FontSize', fontSizeBig);
    c = colorbar;
    set(h,'Interpreter','LaTeX');
    xticks([1:3*day:ncols+1]);
    xticklabels([0:3:num_days+1]);
    xlabel('Time (days)', 'FontSize', fontSizeBig);
    yticks([1:100:num_rows])
    yticklabels([heights(1:100:end)])
    ylim([0, 250])
    xlim([0, train_stop_day])
    ylabel('Height (km)', 'FontSize', fontSizeBig);
    title(sprintf("Wavelet scale %d", ii), 'FontSize', fontSizeBig);
end
hold off;
fig2 = gcf;
fig2.Units = 'inches';
set(fig2, 'position', [0 0 47.7778 18.5555])
if saveFigs
    exportgraphics(fig2, './results/Boulder/fig_wavelet_scale_reconstructions2_2013.eps', ...
        'Resolution', figDPI)
    close(fig2);
end


%% Correlation matrix and adjacency graph
corrs = model.corr_vals;
graph = model.graph;
[nrows, ncols] = size(corrs);
figure
tiledlayout(1, 2, 'TileSpacing','Compact', 'Padding', 'Compact');
nexttile;
hold on
surf(log10(corrs));
shading flat
view([0, 90]);
h = set(gca,'FontSize', fontSize);
colormap('winter')
c = colorbar;
c.FontSize = fontSize;
set(h,'Interpreter','LaTeX');
xlim([1,11])
ylim([1,11])
xticks([1:ncols]);
xlabel('Wavelet index', 'FontSize', fontSize);
yticks([1:nrows]);
ylabel('Wavelet index', 'FontSize', fontSize);
title('Correlation matrix', 'FontSize', fontSize);
hold off;
% fig = gcf;
% fig.Units = 'inches';
% set(fig, 'position', [0 6.1250 14 7.9861])
% if saveFigs
%     exportgraphics(fig, './results/Boulder/fig_adjacency_2013.eps', ...
%         'Resolution', figDPI)
%     close(fig);
% end

% figure
% tiledlayout(1, 1, 'TileSpacing','Compact', 'Padding', 'Compact');
nexttile;
hold on;
plot(graph, 'LineWidth', 4, 'MarkerSize', 8, 'NodeFontSize', fontSize+6, ...
    'EdgeAlpha', 0.5)
h = set(gca,'FontSize', fontSize);
title('Correlation Graph', 'FontSize', fontSize);
hold off;
axis off;
fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 9.1806 20.8472 9.3750])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_correlation_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end


%% Wavelet group reconstructions
lefts = [1,3];
bottoms = [2,3];
wave_groups = model.comps;
[nrows, ncols, ngroups] = size(wave_groups);
tiledlayout(ceil(ngroups/2), 2, 'TileSpacing','Compact', 'Padding', 'Compact');
hold on
for ii=1:ngroups
    nexttile;
    surf(wave_groups(:, :, ii));
    shading flat
    view([0, 90]);
    c = colorbar;
    h = set(gca,'FontSize', fontSizeBig);
%     if sum(ii==lefts)>0
        yticks([1:100:num_rows])
        yticklabels([heights(1:100:end)])
        ylim([0, 250])
        ylabel('Height (km)', 'FontSize', fontSizeBig);
%     else
%         yticks([]);
%         yticklabels([]);
%     end
%     if sum(ii==bottoms)>0
        c.Label.String = 'Freq. (MHz)';
        c.FontSize = fontSizeBig;
        xticks(1:5*day:ncols+1);
        xticklabels([0:5:num_days+1]);
        xlabel('Time (days)', 'FontSize', fontSizeBig);
%     else
%         xticks([]);
%         xticklabels([]);
%     end
    ylim([0, 250])
    xlim([0, train_stop_day])
    title(sprintf("Wavelet Scale %d", ii), 'FontSize', fontSizeBig);
end
hold off;
fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 0 47.7778 18.5555])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_wavelet_group_reconstructions_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end


%% Wavelet group 24-hour averages
wave_group_avgs = model.avgs;
[nrows, ncols, ngroups] = size(wave_group_avgs);
tiledlayout(ceil(ngroups/2), 2, 'TileSpacing','Compact', 'Padding', 'Compact');
hold on
for ii=1:ngroups
    nexttile;
    surf(wave_group_avgs(:, 1:day, ii));
    shading flat
    view([0, 90]);
    c = colorbar;
    h = set(gca,'FontSize', fontSizeBig);
    set(h,'Interpreter','LaTeX');
    xticks([1:8:day+1]);
    xticklabels([0:2:24]);
    xlim([0, day])
    xlabel('Time (hours)', 'FontSize', fontSizeBig);
    yticks([1:100:num_rows])
    yticklabels([heights(1:100:end)])
    ylim([0,num_rows])
    ylabel('Height (km)', 'FontSize', fontSizeBig);
    title(sprintf("Wavelet group %d", ii), 'FontSize', fontSizeBig);
end
hold off;
fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 0 47.7778 18.5555])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_wavelet_group_24-hour_averages_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end


%% Wavelet group 24-hour standard deviations
wave_group_flucs = model.flucs;
[nrows, ncols, ngroups] = size(wave_group_flucs);
tiledlayout(ceil(ngroups/2), 2, 'TileSpacing','Compact', 'Padding', 'Compact');
hold on
for ii=1:ngroups
    nexttile;
    surf(wave_group_flucs(:, 1:day, ii));
    shading flat
    view([0, 90]);
    c = colorbar;
    h = set(gca,'FontSize', fontSizeBig);
    set(h,'Interpreter','LaTeX');
    xticks([1:8:day+1]);
    xticklabels([0:2:24]);
    xlim([0, day])
    xlabel('Time (hours)', 'FontSize', fontSizeBig);
    yticks([1:100:num_rows])
    yticklabels([heights(1:100:end)])
    ylim([0, 250])
    ylabel('Height (km)', 'FontSize', fontSizeBig);
    title(sprintf("Wavelet group %d", ii), 'FontSize', fontSizeBig);
end
hold off;
fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 0 47.7778 18.5555])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_wavelet_group_24-hour_stds_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end


%% Total reconstructions
figure
tiledlayout(2, 1, 'TileSpacing','Compact', 'Padding', 'Compact');
nexttile;
surf(ne);
hold on;
shading flat;
view([0, 90]);
c = colorbar;
c.Label.String = 'Plasma Frequency (MHz)';
caxis([0 11]);
xline(num_days_train*day, 'k--', 'linewidth', 4);
title('Measurement - Dataset 1', 'FontSize', fontSize);
h = set(gca,'FontSize', fontSize);
set(h,'Interpreter','LaTeX')
xticks([1:day:num_cols])
xticklabels([])
xlim([0, num_cols])
ylim([0, num_rows])
xlabel([])
yticks([1:100:num_rows])
yticklabels([heights(1:100:end)])
ylabel('Height (km)', 'FontSize', fontSize)
hold off

nexttile;
surf(model.full_recon);
hold on;
shading flat;
view([0, 90]);
xline(num_days_train*day, 'k--', 'linewidth', 4);
c = colorbar;
c.Label.String = 'Plasma Frequency (MHz)';
caxis([0 11]);
title('Modeled', 'FontSize', fontSize);
h = set(gca,'FontSize', fontSize);
set(h,'Interpreter','LaTeX')
xticks([1:day:num_cols])
xticklabels([0:num_days])
xlim([0, 2400])
ylim([0, 250])
xlabel('Time (days)', 'FontSize', fontSize)
yticks([1:100:num_rows])
yticklabels([heights(1:100:end)])
ylabel('Height (km)', 'FontSize', fontSize)
hold off;

fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 6.1667 22.5417 12.3889])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_recon_with_pred_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end


%% Reconstruction through slices
[nrows, ncols] = size(ne);
figure
h_slice = [200, 150, 100];
slices = [100, 50, 1];
% h_slice = [340, 260, 180, 100];
% slices = [241, 161, 81, 1];
nslices = length(slices);
tiledlayout(nslices, 1, 'TileSpacing','Compact', 'Padding', 'Compact');
for ii=1:nslices
    nexttile;
    model_avg = model.full_recon(slices(ii), :);
    N = size(model_avg, 2);
    x_axis = 1:N;
    plot(ne(slices(ii), :), 'r.', 'markersize', 14);
    hold on;
    plot(x_axis, model_avg, 'b-', 'linewidth', 4);
    xline(num_days_train*day, 'k-', 'linewidth', 4);
    legend('Measured', 'Modeled', 'End train data', 'Location', 'northeast');
    xlim([0, day*(num_days+3.45)]);
    h = set(gca,'FontSize', fontSize);
%     set(h,'Interpreter','LaTeX');
    if ii==nslices
        xticks([1:day:ncols+1]);
        xticklabels([0:1:num_days]);
        xlabel('Time (days)', 'FontSize', fontSize);
    else
        xticks([]);
        xticklabels([]);
    end
    ylabel('Freq. (MHz)', 'FontSize', fontSize);
    title(sprintf("%d km", h_slice(ii)), 'FontSize', fontSize);
end
hold off;
fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 0 30.3750 18.5555])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_avg_model_slices1_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end

figure
h_slice = [350, 300, 250];
slices = [250, 200, 150];
nslices = length(slices);
tiledlayout(nslices, 1, 'TileSpacing','Compact', 'Padding', 'Compact');
for ii=1:nslices
    nexttile;
    model_avg = model.full_recon(slices(ii), :);
    N = size(model_avg, 2);
    x_axis = 1:N;
    plot(ne(slices(ii), :), 'r.', 'markersize', 14);
    hold on;
    plot(x_axis, model_avg, 'b-', 'linewidth', 4);
    xline(num_days_train*day, 'k-', 'linewidth', 4);
    legend('Measured', 'Modeled', 'End train data', 'Location', 'northeast');
    xlim([0, day*(num_days+3.45)]);
    h = set(gca,'FontSize', fontSize);
    set(h,'Interpreter','LaTeX');
    xticks([]);
    xticklabels([]);
    ylabel('Freq. (MHz)', 'FontSize', fontSize);
    title(sprintf("%d km", h_slice(ii)), 'FontSize', fontSize);
end
hold off;
fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 0 30.3750 18.5555])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_avg_model_slices2_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end

%% Large error source
lw = 3;
figure
tiledlayout(1, 1, 'TileSpacing','Compact', 'Padding', 'Compact');
nexttile;
h1 = plot(ne(50, train_stop_day:end), 'ro');
h1.LineWidth = lw;
hold on;
h2 = plot(model.full_recon(50, train_stop_day:end), 'b-');
h2.LineWidth = lw;
h3 = plot(iri.ne(72, train_stop_day:end), 'k-^');
h3.LineWidth = lw;
xticks([0:day:5*day])
xticklabels([20:num_days])
ylim([-1, 9])
xlabel('Time (days)', 'FontSize', fontSize)
ylabel('Plasma Frequency (MHz)', 'FontSize', fontSize)
legend('Data', 'SSDMD', 'IRI')
hold off;
h = set(gca,'FontSize', fontSize);
fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 8.9444 21.6944 9.6112])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_error_phase_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end


%%
% Time dependent error, avg. in height (SSDMD vs IRI vs Daily mean)
lw = 2;
figure
tiledlayout(1, 1, 'TileSpacing','Compact', 'Padding', 'Compact');
nexttile;
plot(mean(abs_per_err_ssdmd, 1), 'b-', 'LineWidth', lw); hold on
plot(mean(abs_per_err_iri, 1), 'r-*', 'LineWidth', 2)
xline(num_days_train*day, 'k--', 'linewidth', 4);
legend('SSDMD', 'IRI')
xticks([1:day:num_cols])
xticklabels([0:num_days])
xlim([0, num_cols])
xlabel('Time (days)', 'FontSize', fontSize)
ylabel('MAPE', 'FontSize', fontSize)
h = set(gca,'FontSize', fontSize);
hold off;
fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 9 25.5556 9.5556])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_error_time_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end


%%
% Height dependent error, avg. in time (SSDMD vs IRI vs Daily median)
lw = 3;
figure
tiledlayout(1, 2, 'TileSpacing','Compact', 'Padding', 'Compact');
nexttile;
plot(mean(abs_per_err_ssdmd(:, 1:train_stop_day), 2), heights, 'b-', 'LineWidth', lw);
hold on
plot(mean(abs_per_err_iri(:, 1:train_stop_day), 2), heights, 'r--', 'LineWidth', lw);
plot(mean(abs_per_err_median(:, 1:train_stop_day), 2), heights, 'g-.', 'LineWidth', lw);
legend('SSDMD', 'IRI', 'Daily Median', 'Location', 'northeast')
xlabel('MAPE', 'FontSize', fontSize)
ylabel('Height (km)', 'FontSize', fontSize)
title('Train')
h = set(gca,'FontSize', fontSize);
hold off;
% fig = gcf;
% fig.Units = 'inches';
% set(fig, 'position', [0 6.5417 10.5000 12.0139])
% if saveFigs
%     exportgraphics(fig, './results/Boulder/fig_error_height_train_2013.eps', ...
%         'Resolution', figDPI)
%     close(fig);
% end
% Height dependent error, avg. in time (SSDMD vs IRI vs Daily median)
% lw = 3;
% figure
% tiledlayout(1, 1, 'TileSpacing','Compact', 'Padding', 'Compact');
nexttile;
plot(mean(abs_per_err_ssdmd(:, train_stop_day:end), 2), heights, 'b-', 'LineWidth', lw);
hold on
plot(mean(abs_per_err_iri(:,train_stop_day:end), 2), heights, 'r--', 'LineWidth', lw);
plot(mean(abs_per_err_median(:,train_stop_day:end), 2), heights, 'g-.', 'LineWidth', lw);
legend('SSDMD', 'IRI', 'Daily Median', 'Location', 'northeast')
xlabel('MAPE', 'FontSize', fontSize)
ylabel('Height (km)', 'FontSize', fontSize)
title('Test')
h = set(gca,'FontSize', fontSize);
hold off;
fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 7.8611 20.4722 10.6945])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_error_height_test_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end


%% Error by hour of day
%%%% Train
mape_hours_ssdmd1 = zeros(num_rows, day);
for ii=1:day
    mape_hours_ssdmd1(:, ii) = mean(abs_per_err_ssdmd(:, ii:day:train_stop_day), 2);
end
mape_hours_ssdmd1 = circshift(mape_hours_ssdmd1, 24, 2);
figure
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
surf(mape_hours_ssdmd1);
shading flat
view([0, 90])
% colormap(flipud(bone))
colormap('jet')
cb = colorbar;
cb.Label.String = 'MAPE';
cb.FontSize = fontSize;
caxis([0, 100])
xticks([1:8:day])
xticklabels([mod(([0:2:24]), 24)])
yticks([1:50:250])
yticklabels([100:50:350])
xlim([0, day])
xlabel('Hour of Day (Local)', 'FontSize', fontSize)
ylabel('Height (km)', 'FontSize', fontSize)
title('Train')
h = set(gca,'FontSize', fontSize);
hold off;
% fig = gcf;
% fig.Units = 'inches';
% set(fig, 'position', [-7.9722 17.3194 21.7917 13])
% if saveFigs
%     exportgraphics(fig, './results/Boulder/fig_error_hours_train_2013.eps', ...
%         'Resolution', figDPI)
%     close(fig);
% end

%%%% Test
mape_hours_ssdmd2 = zeros(num_rows, day);
for ii=1:day
    mape_hours_ssdmd2(:, ii) = mean(abs_per_err_ssdmd(:, ii+train_stop_day:day:end), 2);
end
mape_hours_ssdmd2 = circshift(mape_hours_ssdmd2, 24, 2);
% figure
% tiledlayout(1, 1, 'TileSpacing','Compact', 'Padding', 'Compact');
nexttile;
surf(mape_hours_ssdmd2);
shading flat
view([0, 90])
colormap('jet')
cb = colorbar;
cb.Label.String = 'MAPE';
cb.FontSize = 36;
caxis([0, 100])
xticks([1:8:day])
xticklabels([mod(([0:2:24]), 24)])
yticks([1:50:250])
yticklabels([100:50:350])
xlim([0, day])
xlabel('Hour of Day (Local)', 'FontSize', 36)
ylabel('Height (km)', 'FontSize', 36)
title('Test')
h = set(gca,'FontSize', 36);
hold off;
fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 7.3472 36.3333 11.2083])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_error_hours_test_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end


%% Error histogram
norm_type = 'pdf';

%%%%% Fit
pt_err_ssdmd_fit = ne(:, 1:train_stop_day) - ...
    model.full_recon(:, 1:train_stop_day);
pt_err_iri_fit = ne(:, 1:train_stop_day) - ...
    iri.ne(:, 1:train_stop_day);

figure
tiledlayout(1, 2, 'TileSpacing','Compact', 'Padding', 'Compact');
nexttile;
h2 = histogram(pt_err_iri_fit);
h2.Normalization = norm_type;
h2.FaceColor = 'r';
h2.FaceAlpha = 0.5;
h2.EdgeAlpha = 0.5;
hold on
h1 = histogram(pt_err_ssdmd_fit);
h1.Normalization = norm_type;
h1.FaceColor = 'b';
h1.FaceAlpha = 0.5;
h1.EdgeAlpha = 0.5;
hold off
legend('IRI', 'SSDMD')
xlabel('$${\bf Y} - {\bf \hat{Y}}$$ (MHz)', 'Interpreter', 'latex')
ylabel('PDF')
title('Train')
h = set(gca,'FontSize', fontSize);
% fig = gcf;
% fig.Units = 'inches';
% set(fig, 'position', [0 10.5556 11.6944 8.0000])
% if saveFigs
%     exportgraphics(fig, './results/Boulder/fig_error_hist_train_2013.eps', ...
%         'Resolution', figDPI)
%     close(fig);
% end

%%%%% Test
pt_err_ssdmd_test = ne(:, train_stop_day:end) - ...
    model.full_recon(:, train_stop_day:end);
pt_err_iri_test = ne(:, train_stop_day:end) - ...
    iri.ne(:, train_stop_day:end);

% figure
% tiledlayout(1, 1, 'TileSpacing','Compact', 'Padding', 'Compact');
nexttile;
h2 = histogram(pt_err_iri_test);
h2.Normalization = norm_type;
h2.FaceColor = 'r';
h2.FaceAlpha = 0.5;
h2.EdgeAlpha = 0.5;
hold on
h1 = histogram(pt_err_ssdmd_test);
h1.Normalization = norm_type;
h1.FaceColor = 'b';
h1.FaceAlpha = 0.5;
h1.EdgeAlpha = 0.5;
hold off
legend('IRI', 'SSDMD')
xlabel('$${\bf Y} - {\bf \hat{Y}}$$ (MHz)', 'Interpreter', 'latex')
ylabel('PDF')
title('Test')
h = set(gca,'FontSize', fontSize);
fig = gcf;
fig.Units = 'inches';
set(fig, 'position', [0 9.5972 22.4722 8.9584])
if saveFigs
    exportgraphics(fig, './results/Boulder/fig_error_hist_test_2013.eps', ...
        'Resolution', figDPI)
    close(fig);
end


