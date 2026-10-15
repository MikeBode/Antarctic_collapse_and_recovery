function Figure2_ModelledObserved()
% FIGURE2_MODELLEDOBSERVED  Draws Figure 2.
%
% Compares the best-fitting model prediction with the observed annual
% proportion of Measures entering into force, and with the moving average
% from Figure 1. Needs FitModel and Figure1_MeasuresAndNetworks to have run.
%
% Output: Figures/Figure2_ModelledObserved.tif

F  = load(fullfile('CreatedData','BestFit.mat'), 'years', 'pkBestGlobal', 'years_obs', 's_u', 'n_u');
F1 = load(fullfile('CreatedData','Figure1_MovingAverage.mat'), 'xxs', 'yys');
years = F.years;

% Figure and colours (matched to Figure 1)
figure(2); clf;
subplot('position', [0.07 0.09 0.93 0.88]); hold on
FS = 16;
set(gca, 'fontsize', FS-2); box on
Col     = slanCM('ice', 20);
LineCol = Col(17,:);   % observed data
gridCol = Col(18,:);
FontCol = Col(7,:);    % model line and text

% Leave out the first two and last two years
plot_years = 3:length(years)-2;

% Moving average from Figure 1
plot(F1.xxs, F1.yys, '-', 'Color', LineCol, 'LineWidth', 3);

% Observed annual proportions
for p = 1:length(plot_years)
    k = plot_years(p);
    plot(F.years_obs(k), F.s_u(k) ./ F.n_u(k), '.', 'MarkerSize', 12, 'Color', LineCol);
end

% Best-fit model prediction, smoothed for display
x_model = years(plot_years);
y_model = F.pkBestGlobal(plot_years);
xx = linspace(min(x_model), max(x_model), 400);
yy = interp1(x_model, y_model, xx, 'pchip');
plot(xx, yy, '-', 'Color', FontCol, 'LineWidth', 2.25);
plot(x_model, y_model, '.', 'Markersize', 11, 'Color', FontCol);

% Legend
TY = 2006;
Th = 0.25;
plot(TY+[0 4], [1 1].*Th, '-', 'Color', FontCol, 'LineWidth', 3);
plot(TY+2, Th, '.', 'Color', FontCol, 'Markersize', 14);
text(TY+4.5, Th, 'FJ Model prediction', 'fontsize', FS-4, 'Color', FontCol);
Th = Th - 0.075;
plot(TY+[0 4], [1 1].*Th, '-', 'Color', LineCol, 'LineWidth', 3);
text(TY+4.5, Th, 'ATCM obs (moving average)', 'fontsize', FS-4, 'Color', FontCol);
Th = Th - 0.075;
plot(TY+2, Th, '.', 'MarkerSize', 18, 'Color', LineCol);
text(TY+4.5, Th, 'ATCM obs (annual probability)', 'fontsize', FS-4, 'Color', FontCol);

% Hand-drawn axes and grid
ylim([0 1.025]);
xlim([min(years)+1, max(years)+1]);
xB = [1960 2025];
for yr = 1970:10:2020
    text(yr, -0.075, num2str(yr), 'horizontalalignment','center', 'fontsize', FS-4, 'color', FontCol)
end
for gy = 0:0.2:1
    if gy == 0.2
        plot([1960 2005], [gy gy], '--', 'linewidth', 0.5, 'color', gridCol)   % stops short of the legend
    else
        plot(xB, [gy gy], '--', 'linewidth', 0.5, 'color', gridCol)
    end
    text(xB(1)+1.5, gy, sprintf('%.0f', gy*100), 'horizontalalignment','right', ...
        'verticalalignment','middle', 'fontsize', FS-4, 'color', FontCol)
end
text(xB(1)-2, 0.5, 'Probability a Measure enters into force', 'rotation', 90, ...
    'horizontalalignment','center', 'fontsize', FS-2, 'color', FontCol)
axis off

% Export
DIM = [0 0 35 14]*0.75;
MakeTIFF(fullfile('Figures','Figure2_ModelledObserved'), DIM);
fprintf('  Saved Figures/Figure2_ModelledObserved.tif\n');
end
