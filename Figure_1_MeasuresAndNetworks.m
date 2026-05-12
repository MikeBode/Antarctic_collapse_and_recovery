%% Figure_S1_DataAndNetworks_likeOriginal.m
% Recreate the “network + consensus probability” figure in the original style,
% but using your newer data sources and with an option to use span-averaged
% networks or raw per-year networks.
%
% - Top panel: small “network glyphs” for selected ATCMs laid out across time
% - Bottom panel: observed probability of reaching consensus per year,
%                 shown as a smoothed line and high-likelihood dots from Beta posteriors
%
% Data inputs:
%   * Span-averaged networks (if useSpan = true):
%       CreatedData/AveragedNetworks_[No|Re]visions/Collaboration_by_year_span_S1_S2.mat
%   * Raw networks (if useSpan = false):
%       Data/Collaboration_by_year.mat
%   * ATCM measures (counts):
%       CreatedData/MeasureData_NoSATCM.mat   (expects MeasuresInfo)
%   * Optional Roman numerals for labels:
%       ../OLD/V1/Data/ATCM_roman_numerals.xlsx
%
% Output:
%   * On-screen figure (dark blue theme like original)
%   * Saves Likelihood and Observed trend to Data/Collapse_trend_[RevFlag].mat
%   * Optional: TIFF export at the bottom (commented)

clear; clc;

% ------------------- USER CHOICES -------------------
Span    = [3 2];         % [S1 S2] used only if useSpan=true
Span    = [0 0];         % [S1 S2] used only if useSpan=true
showIdx = [3 7 11 16 22 30 39];   % which ATCMs to show as network glyphs along the top
smoothSpan = 7;          % LOWESS span (odd integer recommended) for the bottom trend
dotThresh  = 0.75;       % probability density threshold for drawing beta dots
% ----------------------------------------------------

% Load measure counts (expects MeasuresInfo: [year; n_trials; n_success])
load('CreatedData/MeasureData_NoSATCM','MeasuresInfo');

% Roman numerals for labels (optional)
[~, t_dates] = xlsread('../OLD/V1/Data/ATCM_roman_numerals.xlsx');

PlotVer = 1;
%% Version 1: smoothed line + dots
if PlotVer == 1
    for IncludeRevisions = 1
        % ---------- Load networks ----------
        if IncludeRevisions == 0
            fnet = sprintf('CreatedData/AveragedNetworks_NoRevisions/Collaboration_by_year_span_%d_%d.mat', Span(1), Span(2));
        else
            fnet = sprintf('CreatedData/AveragedNetworks_Revisions/Collaboration_by_year_span_%d_%d.mat', Span(1), Span(2));
        end
        S = load(fnet,'WP_Collaboration');
        WP_Collaboration = S.WP_Collaboration;
        Years    = cell2mat(WP_Collaboration(:,1));

        % ---------- FIGURE + STYLE ----------
        figure(1); clf; FS = 18;

        % === NETWORK SUBPLOT (top strip) ===
        axNet = subplot('position',[0.07 0.72 0.91 0.25]); hold(axNet,'on');
        xlim(axNet, [min(Years)-1 max(Years)+1]);
        ylim(axNet, [1.05 1.5]);
        axis(axNet,'off');

        % Colour scheme
        Col = slanCM('ice',20);
        NetworkCL = Col([17 14 9 5],:);
        LineCol = Col(17,:);
        MarkCol = Col(9,:);
        gridCol = Col(18,:);
        FontCol = Col(7,:);

        % Parameters for the mini network glyphs (kept very close to your original)
        r       = 2.5.*[1 0.05];  % glyph x/y scaling
        YS      = 1.3;            % y at which glyphs are centered
        MKT     = Col(7,:);           % node marker
        MKS     = 10;              % marker size
        LW      = [1 1];      % [thin, thick] edge widths
        FSg     = 10;             % glyph label font size
        TShift  = 0.175;          % vertical offset for glyph label

        % Clamp indices to valid range
        showIdx = showIdx(showIdx>=1 & showIdx<=size(WP_Collaboration,1));

        % Draw selected ATCM networks spaced by their actual year on the x-axis
        for y = showIdx(:).'
            ThisYear = WP_Collaboration{y,1};
            TitleStr = t_dates{y,1};
            Q = add_network_to_figure(WP_Collaboration{y,4}, ...  % adjacency
                ThisYear, r, YS, TitleStr, ...             % x, glyph size, y, label text
                TShift, FSg, LW, MKT, MKS, ...             % label shift/size, widths, node style
                NetworkCL, 0);                                        % W=0 (transparent badge)
        end
        YLB = [2021 2023]+1;
        s = 1.2+1.45.*[0 0.05 0.1 0.15];
        for i = 1:4
            plot(YLB,[0 0]+s(i),'-','linewidth',1.5,'color',NetworkCL(i,:))
            text(YLB(2)+0.2,s(i),sprintf('$%d\\%%$', round(100*Q(i))),'fontsize',FS-6,'verticalalignment','middle')
        end

        % === PROBABILITY SUBPLOT (bottom) ===
        axProb = subplot('position',[0.07 0.07 0.91 0.60]); hold(axProb,'on');
        ylabel(axProb,{'Proportion of Measures in force'},'FontSize',FS-1);

        xB = [1958 2025];
        for y = [1960:10:2020]
            text(y,-0.075,num2str(y),...
                'horizontalalignment','center',...
                'fontsize',FS-4,'color',FontCol)
            % plot([y y],[0 1],'--','linewidth',0.5,'color',gridCol)
        end
        for yy = 0:0.2:1
            plot(xB,[yy yy],'--','linewidth',0.5,'color',gridCol)
            % text(xB(1)-0.5,yy,sprintf('%.1f', yy),'horizontalalignment','right','verticalalignment','middle',...
            %     'fontsize',FS-4,'color',FontCol)
            text(xB(1)-0.5,yy,sprintf(['%.0f'], yy*100),'horizontalalignment','right','verticalalignment','middle',...
                'fontsize',FS-4,'color',FontCol)
        end
        text(xB(1)-4,0.5,'Percent of Measures in force','rotation',90',...
            'horizontalalignment','center',...
            'fontsize',FS-2,'color',FontCol)

        text(xB(1)-6,1.05,'(B)','fontsize',FS-2,'color',FontCol)
        text(xB(1)-6,1.45,'(A)','fontsize',FS-2,'color',FontCol)
        axis off

        % ---------- Build per-year beta likelihood & expected proportion ----------
        years_m   = MeasuresInfo(1,:);          % year
        n_trials  = MeasuresInfo(2,:);          % trials
        n_success = MeasuresInfo(3,:);          % successes

        probGrid  = linspace(0,1,100);          % grid for the beta pdf (fine)
        probGridD = linspace(0,1,30);           % coarser grid for plotting dots

        Likelihood = cell(numel(years_m),5);
        expect     = zeros(1,numel(years_m));

        for k = 1:numel(years_m)-1
            % Beta posterior parameters (flat prior -> +1, +1)
            a = n_success(k) + 1;
            b = (n_trials(k) - n_success(k)) + 1;

            pdf = betapdf(probGrid, a, b);
            pdf = pdf ./ sum(pdf);  % normalized for visualization (optional)

            expect(k) = n_success(k) / max(n_trials(k),1);  % observed proportion
            plot(years_m(k), expect(k), '.', 'color', MarkCol,'MarkerSize', 15);

        end

        % ---------- Smooth the observed proportion (LOWESS), like your original ----------
        prop_m = smooth(expect, 4);

        % ----- Curve with pchip interpolation between -----
        x = years_m(3:end-3);
        y = prop_m(3:end-3);
        xx = linspace(min(x), max(x), 400);            % finer grid for smoothness
        yy = interp1(x, y, xx, 'pchip');               % full pchip interpolation
        hLine = plot(xx, yy, '-', 'Color', LineCol, 'LineWidth', 2.5);
        uistack(hLine, 'bottom');
        set(axProb,'FontSize',FS-2,'layer','top');
        xlabel(axProb,'Year','FontSize',FS);
        ylim([0 1]), xlim([1957 2025])

        TY = 2008; Th = 0.1;
        plot(TY+[0 4],[1 1].*Th,'-', 'Color', LineCol, 'LineWidth', 4);
        text(TY+4.5,Th,'Moving average','fontsize',FS, 'Color', LineCol.^4);

        DIM = [0 0 35 17];
        if IncludeRevisions == 0
            Make_TIFF(['Figures/Figure1_MeasuresAndNetworks_NoRevisions'], DIM);
        else
            Make_TIFF(['Figures/Figure1_MeasuresAndNetworks_Revisions'], DIM);
        end
    end
end


%% Version 2: bar + dots
if PlotVer == 2
fnet = sprintf('CreatedData/AveragedNetworks_Revisions/Collaboration_by_year_span_%d_%d.mat', Span(1), Span(2));
S = load(fnet,'WP_Collaboration');
WP_Collaboration = S.WP_Collaboration;
Years    = cell2mat(WP_Collaboration(:,1));

% ---------- FIGURE + STYLE ----------
figure(1); clf; FS = 16;

% === NETWORK SUBPLOT (top strip) ===
axNet = subplot('position',[0.07 0.72 0.91 0.25]); hold(axNet,'on');
xlim(axNet, [min(Years)-1 max(Years)+1]);
ylim(axNet, [1.05 1.5]);
axis(axNet,'off');

% Colour scheme
Col = slanCM('ice',20);
NetworkCL = Col([17 14 9 5],:);
LineCol = Col(17,:);
MarkCol = Col(9,:);
gridCol = Col(18,:);
FontCol = Col(7,:);

% Parameters for the mini network glyphs (kept very close to your original)
r       = 2.5.*[1 0.05];  % glyph x/y scaling
YS      = 1.3;            % y at which glyphs are centered
MKT     = Col(7,:);           % node marker
MKS     = 10;              % marker size
LW      = [1 1];      % [thin, thick] edge widths
FSg     = 10;             % glyph label font size
TShift  = 0.175;          % vertical offset for glyph label

% Clamp indices to valid range
showIdx = showIdx(showIdx>=1 & showIdx<=size(WP_Collaboration,1));

% Draw selected ATCM networks spaced by their actual year on the x-axis
for y = showIdx(:).'
    ThisYear = WP_Collaboration{y,1};
    TitleStr = t_dates{y,1};
    Q = add_network_to_figure(WP_Collaboration{y,4}, ...  % adjacency
        ThisYear, r, YS, TitleStr, ...             % x, glyph size, y, label text
        TShift, FSg, LW, MKT, MKS, ...             % label shift/size, widths, node style
        NetworkCL, 0);                                        % W=0 (transparent badge)
end
YLB = [2021 2023]+1;
s = 1.2+1.45.*[0 0.05 0.1 0.15];
for i = 1:4
    plot(YLB,[0 0]+s(i),'-','linewidth',1.5,'color',NetworkCL(i,:))
    text(YLB(2)+0.2,s(i),sprintf('$%d\\%%$', round(100*Q(i))),'fontsize',FS-6,'verticalalignment','middle')
end

% === PROBABILITY SUBPLOT (bottom) ===
axProb = subplot('position',[0.07 0.07 0.91 0.60]); hold(axProb,'on');
ylabel(axProb,{'Proportion of Measures in force'},'FontSize',FS-1);

xB = [1958 2025];
for y = [1960:10:2020]
    text(y,-0.075,num2str(y),...
        'horizontalalignment','center',...
        'fontsize',FS-4,'color',FontCol)
end
for yy = 0:0.2:1
    plot(xB,[yy yy],'--','linewidth',0.5,'color',gridCol)
    text(xB(1)-0.5,yy,sprintf(['%.0f'], yy*100),'horizontalalignment','right','verticalalignment','middle',...
        'fontsize',FS-4,'color',FontCol)
end
text(xB(1)-4,0.5,'Percent of Measures in force','rotation',90',...
    'horizontalalignment','center',...
    'fontsize',FS-2,'color',FontCol)

text(xB(1)-6,1.05,'(B)','fontsize',FS+2,'color',FontCol)
text(xB(1)-6,1.45,'(A)','fontsize',FS+2,'color',FontCol)
axis off

% ---------- Build per-year beta likelihood & expected proportion ----------
years_m   = MeasuresInfo(1,:);          % year
n_trials  = MeasuresInfo(2,:);          % trials
n_success = MeasuresInfo(3,:);          % successes

probGrid  = linspace(0,1,100);          % grid for the beta pdf (fine)
probGridD = linspace(0,1,30);           % coarser grid for plotting dots

Likelihood = cell(numel(years_m),5);
expect     = zeros(1,numel(years_m));

for k = 1:numel(years_m)-1
    expect(k) = n_success(k) / max(n_trials(k),1);  % observed proportion
end

for k = 1:numel(years_m)-1
    if k > 1
        pp = patch(years_m(k)+[0 0 1 1].*(years_m(k+1)-years_m(k)),...
            [0 expect(k) expect(k) 0],LineCol);
        set(pp,'edgecolor','none','facecolor',LineCol)
    
        plot([years_m(k-1) years_m(k)],[expect(k-1) expect(k)],'k')
    elseif k == 1
        pp = patch(years_m(k)+[0 0 1 1],...
            [0 expect(k) expect(k) 0],LineCol);
        set(pp,'edgecolor','none','facecolor',LineCol)
    end
end

% % ---------- Smooth the observed proportion (LOWESS), like your original ----------
% prop_m = smooth(expect, 4);
% hLine = plot(years_m(3:end-3), prop_m(3:end-3), ...
%     '-', 'Color', LineCol, 'LineWidth', 3);
% uistack(hLine, 'bottom');

set(axProb,'FontSize',FS-2,'layer','top');
xlabel(axProb,'Year','FontSize',FS);
ylim([0 1]), xlim([1957 2025])

DIM = [0 0 35 17];
Make_TIFF(['Figures/Figure1_MeasuresAndNetworks_Revisions_ver2'], DIM);
end
