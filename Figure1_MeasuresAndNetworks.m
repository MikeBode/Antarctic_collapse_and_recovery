function Figure1_MeasuresAndNetworks(Span)
% FIGURE1_MEASURESANDNETWORKS  Draws Figure 1.
%
% (A) The co-authorship network for every ATCM, laid out in rows that
%     zigzag forward in time.
% (B) The observed proportion of Measures entering into force at each ATCM,
%     with a smoothed moving average.
%
% Input
%   Span  [S1 S2] span of the networks to draw
%
% Outputs
%   Figures/Figure1_MeasuresAndNetworks.tif
%   CreatedData/Figure1_MovingAverage.mat (the moving average line, used by Figure 2)

load(fullfile('CreatedData','MeasureData.mat'), 'MeasuresInfo');
[d_dates, ~] = xlsread(fullfile('RawData','ATCM_roman_numerals.xlsx'));

fnet = fullfile('CreatedData','AveragedNetworks', ...
    sprintf('Collaboration_by_year_span_%d_%d.mat', Span(1), Span(2)));
S = load(fnet, 'WP_Collaboration');
WP_Collaboration = S.WP_Collaboration;
Years    = cell2mat(WP_Collaboration(:,1));
NumATCMs = size(WP_Collaboration, 1);

% The network labels come from the meeting list, so check it lines up
assert(isequal(d_dates(1:NumATCMs,3), Years), ...
    'Years in ATCM_roman_numerals.xlsx do not match the network years.');

% Figure and colours
figure(1); clf; FS = 18;
Col = slanCM('ice', 20);
NetworkCL = Col([17 14 9 5],:);
LineCol = Col(17,:);
MarkCol = Col(9,:);
gridCol = Col(18,:);
FontCol = Col(7,:);

% ---------- Panel A: networks ----------
YX = 0.35;
subplot('position', [0.03 YX-0.02 0.95 0.97-YX]); hold on;
PerRow  = 7;
NumRows = ceil(NumATCMs/PerRow);
SC      = [4 3.5];

% Glyph positions (x, y) and the connecting line (xN, yN). Odd rows run
% left to right, even rows right to left.
x = []; y = [];
xN = []; yN = [];
for i = 1:NumRows
    for k = 1:PerRow
        if mod(i,2) == 0
            x = [x PerRow+1-k];
        else
            x = [x k];
        end
        y = [y i];
    end
    for j = 1:100
        if mod(i,2) == 0
            xN = [xN PerRow-(j-1)/98*PerRow];
        else
            xN = [xN (j-2)/98*PerRow];
        end
        yN = [yN i];
    end
end
x  = (x - 0.5) * SC(1);
xN = xN * SC(1);
y  = y  * SC(2);
yN = yN * SC(2);

% Flip so the earliest row is at the top
yMax = max(y);
y  = ((yMax + 1) - y/SC(2))  * SC(2);
yN = ((yMax + 1) - yN/SC(2)) * SC(2);

% Connecting line behind the networks
plot(xN(15:end-30), yN(15:end-30), '-', 'linewidth', 3, 'color', LineCol)

% Arrow heads where the line turns at the end of each row
as = 0.2; bs = 0.5;
for i = 1:NumRows-1
    yRow = (yMax + 1 - i)*SC(2);
    if mod(i,2) == 1
        plot(PerRow*SC(1)+[-as 0 as], yRow - SC(2)/2 + [0 -bs 0], '-', 'linewidth', 3, 'color', LineCol)
    else
        plot([-as 0 as]-0.3, yRow - SC(2)/2 + [0 -bs 0], '-', 'linewidth', 3, 'color', LineCol)
    end
end

% One network glyph per ATCM, labelled with its year
r      = [1 1]*1.25;   % glyph radius
LW     = [0.75 1];     % line widths for light and heavy edges
MKT    = Col(7,:);     % node and label colour
MKS    = 9;            % node size
FSg    = 14;           % label font size
TShift = 0;            % label offset from glyph centre
for i = 1:NumATCMs
    TitleStr = num2str(d_dates(i,3));
    AddNetworkToFigure(WP_Collaboration{i,4}, x(i), r, y(i), TitleStr, TShift, FSg, LW, MKT, MKS, NetworkCL);
end
drawnow

text(-2.1, (yMax+1)*SC(2)-3, '(A)', 'fontsize', FS-2, 'color', FontCol);
axis off
axis tight
XL = xlim; W = XL(2)-XL(1);
xlim(XL + [-W W]./50)

% ---------- Panel B: proportion in force ----------
axProb = subplot('position', [0.07 0.07 0.91 YX-0.07-0.05]); hold(axProb, 'on');

% Hand-drawn axes and grid
xB = [1958 2025];
for yr = 1960:10:2020
    text(yr, -0.075, num2str(yr), 'horizontalalignment','center', 'fontsize', FS-4, 'color', FontCol)
end
for yy = 0:0.2:1
    plot(xB, [yy yy], '--', 'linewidth', 0.5, 'color', gridCol)
    text(xB(1)-0.5, yy, sprintf('%.0f', yy*100), 'horizontalalignment','right', ...
        'verticalalignment','middle', 'fontsize', FS-4, 'color', FontCol)
end
text(xB(1)-4, 0.5, 'Percent of Measures in force', 'rotation', 90, ...
    'horizontalalignment','center', 'fontsize', FS-2, 'color', FontCol)
text(xB(1)-6.5, 1.05, '(B)', 'fontsize', FS-2, 'color', FontCol)
axis off

% Observed proportion for each year. The final year is left at zero and
% not plotted.
years_m   = MeasuresInfo(1,:);
n_trials  = MeasuresInfo(2,:);
n_success = MeasuresInfo(3,:);
expect = zeros(1, numel(years_m));
for k = 1:numel(years_m)-1
    expect(k) = n_success(k) / max(n_trials(k), 1);
    plot(years_m(k), expect(k), '.', 'color', MarkCol, 'MarkerSize', 15);
end

% Moving average, trimmed at the ends and interpolated for a smooth line
prop_m = smooth(expect, 4);
xs  = years_m(3:end-3);
ys  = prop_m(3:end-3);
xxs = linspace(min(xs), max(xs), 400);
yys = interp1(xs, ys, xxs, 'pchip');
hLine = plot(xxs, yys, '-', 'Color', LineCol, 'LineWidth', 3);
uistack(hLine, 'bottom');

save(fullfile('CreatedData','Figure1_MovingAverage.mat'), 'xxs', 'yys', 'years_m', 'expect');

set(axProb, 'FontSize', FS-2, 'layer', 'top');
ylim(axProb, [0 1.025]); xlim(axProb, [1957 2025]);

% Legend
TY = 2008; Th = 0.14; FSl = FS - 3;
plot(TY+[0 4], [1 1]*Th, '-', 'Color', LineCol, 'LineWidth', 4);
text(TY+4.5, Th, 'Moving average', 'fontsize', FSl, 'Color', LineCol.^4);
Th = Th + 0.125;
plot(TY+2, Th, '.', 'color', MarkCol, 'MarkerSize', 15);
text(TY+4.5, Th, 'ATCM observed', 'fontsize', FSl, 'Color', LineCol.^4);

% Export
DIM = [0 0 30 38].*0.8;
MakeTIFF(fullfile('Figures','Figure1_MeasuresAndNetworks'), DIM);
fprintf('  Saved Figures/Figure1_MeasuresAndNetworks.tif\n');
end
