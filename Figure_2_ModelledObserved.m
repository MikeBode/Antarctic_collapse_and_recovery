% Figure_2_ModelledObserved()
% Bootstrap-fit against a flat matrix of simulated pk time series,
% with OPTIONAL restrictions on the span parameters S1 (past) and S2 (future).
%
% Flat store: Data/pk_flat_store.mat with:
%   Mflat: (4+Y_flat) x N  where rows are [R; B; S1; S2; pk(1:Y)]
%   years: Y_flat x 1      (ATCM years)
%
% Observations: Data/Likelihood.mat (per-year counts -> rates)
%
% Loss metric:   'LL'  (binomial log-likelihood from counts)
%                'SSD' (sum of squared differences vs observed rate)
%
% Output: figure with posterior heat rectangles + bootstrap quantile envelope.

clear;

for m = 2
    if m == 1
        metric = 'LL';
    elseif m == 2
        metric = 'SSD';
    end

    for IncludeRevisions = 1
        %% PREP
        load CreatedData/MeasureData
        if IncludeRevisions == 0
            file_append = 'NoRevisions';
        elseif IncludeRevisions == 1
            file_append = 'Revisions';
        end
        eval(['load CreatedData/Timeseries_' file_append '/pk_flat_store Mflat years'])
        years_flat  = years(:);                   % Y_flat x 1
        Y_flat      = numel(years_flat);

        % ---- Span restrictions (inclusive). Leave [] to disable a bound. ----
        SpanPastRange   = [0 7];                 % limits for S1 (past window). [] => no limit
        SpanFutureRange = [0 7];                 % limits for S2 (future window). [] => no limit

        % Bootstrap
        B = 200;                                   % number of bootstrap resamples

        % Fit window on the aligned year axis
        fitYears  = [1964 2020];                    % [] to use all aligned years
        dropZero  = true;                           % drop years with zero trials

        % Plot appearance
        lineWidth = 2;
        lineColor = [0 0 0];                        % quantile envelope color (black)
        redWash   = [linspace(1,1,256)', linspace(1,0.7,256)', linspace(1,0.7,256)']; % white→soft red

        % Extract parameters
        R_all  = double(Mflat(1,:));
        B_all  = double(Mflat(2,:));
        S1_all = double(Mflat(3,:));
        S2_all = double(Mflat(4,:));

        % Apply span restrictions
        keep = true(1, size(Mflat,2));
        if ~isempty(SpanPastRange)
            keep = keep & (S1_all >= SpanPastRange(1)) & (S1_all <= SpanPastRange(2));
        end
        if ~isempty(SpanFutureRange)
            keep = keep & (S2_all >= SpanFutureRange(1)) & (S2_all <= SpanFutureRange(2));
        end
        if ~any(keep)
            error('No parameter sets remain after span filtering. Adjust SpanPastRange/SpanFutureRange.');
        end

        % Filter flat store
        Mflat   = Mflat(:, keep);
        R_all   = R_all(keep);
        B_all   = B_all(keep);
        S1_all  = S1_all(keep);
        S2_all  = S2_all(keep);
        Pk_all  = double(Mflat(5:end,:));           % Y_flat x Ncols (filtered)
        Ncols   = size(Pk_all,2);

        fprintf('Span-restricted search over %d parameter sets (%.1f%% of original).\n', ...
            Ncols, 100*Ncols/size(Mflat,2));

        %% -------------------- LOAD OBSERVATIONS (COUNTS -> RATE) --------------------
        y_obs_raw = cell2mat(Likelihood(:,1));
        a_raw     = cell2mat(Likelihood(:,4));
        b_raw     = cell2mat(Likelihood(:,5));
        s_raw = a_raw - 1;                    % successes
        n_raw = a_raw + b_raw - 2;            % trials

        % Aggregate by year (handles potential duplicates)
        [years_obs, ~, g] = unique(y_obs_raw(:));
        s_u = accumarray(g, s_raw(:), [], @sum);
        n_u = accumarray(g, n_raw(:), [], @sum);
        p_u = s_u ./ max(n_u, 1);

        %% -------------------- ALIGN YEARS --------------------
        [years, iFlat, iObs] = intersect(years_flat, years_obs, 'stable');
        if isempty(years)
            error('No common years between flat store and observations.');
        end

        Pk   = Pk_all(iFlat, :);              % Y x Ncols
        s    = double(s_u(iObs));             % Y x 1
        n    = double(n_u(iObs));             % Y x 1
        pobs = double(p_u(iObs));             % Y x 1
        Y    = numel(years);

        assert(size(Pk,1)==Y && numel(s)==Y && numel(n)==Y, 'Alignment mismatch.');

        %% -------------------- FIT MASK (BUILT ON ALIGNED YEARS) --------------------
        fitMask = true(Y,1);
        if ~isempty(fitYears)
            fitMask = fitMask & years>=fitYears(1) & years<=fitYears(2);
        end
        if dropZero
            fitMask = fitMask & (n>0);
        end
        assert(any(fitMask), 'No fitting years after applying bounds/filters.');
        Yf = sum(fitMask);

        %% -------------------- PRECOMPUTE LOSS CONTRIBUTIONS --------------------
        switch upper(metric)
            case 'LL'
                Pfit   = min(max(Pk(fitMask,:), 1e-12), 1-1e-12);            % avoid log(0)
                Lmat   = (s(fitMask).*log(Pfit)) + ((n(fitMask)-s(fitMask)).*log1p(-Pfit)); % Yf x Ncols
            case 'SSD'
                Diff2  = (Pk(fitMask,:) - pobs(fitMask)).^2;                  % Yf x Ncols
            otherwise
                error('metric must be ''LL'' or ''SSD''.');
        end

        %% -------------------- LOSS LANDSCAPES: R×B for each Span = [S1 S2] --------------------
        plot_landscape = 0;
        if plot_landscape == 1

            % Collapse per-column loss on the fitting window to a single scalar
            switch upper(metric)
                case 'LL'
                    colScore = sum(Lmat, 1);          % higher is better
                    cbarLabel = 'Sum log-likelihood';
                case 'SSD'
                    colScore = -sum(Diff2, 1);        % negate so that higher is better for plotting
                    cbarLabel = '- Sum of squared differences';
            end
            colScore = double(colScore(:));          % Ncols x 1

            % Restrict to Span ranges shown on panels
            S1_vals = 0:7;
            S2_vals = 0:7;

            % Unique R and B grids present after filtering
            R_vals  = unique(R_all(:)).';
            B_vals  = unique(B_all(:)).';
            nRx = numel(R_vals);
            nBy = numel(B_vals);

            for ii = 1:numel(S1_vals)
                for jj = 1:numel(S2_vals)
                    mask = (S1_all == S1_vals(ii)) & (S2_all == S2_vals(jj));
                    if ~any(mask)
                        continue
                    end
                    figure(1), clf

                    % Map to grid indices
                    [~, ix] = ismember(R_all(mask), R_vals);
                    [~, iy] = ismember(B_all(mask), B_vals);

                    % Ensure column vectors and valid subscripts
                    ix = ix(:);
                    iy = iy(:);
                    vals = colScore(mask);
                    vals = vals(:);
                    ok = ix > 0 & iy > 0;   % guard against any non-matches

                    % Best score per (R,B) cell, NaN for empty cells
                    Z = accumarray([iy(ok) ix(ok)], vals(ok), [nBy nRx], @max, NaN);
                    [RR, BB] = meshgrid(R_vals, B_vals);

                    [bestPanelScore, linBest] = max(Z(:));
                    [iyBest, ixBest] = ind2sub(size(Z), linBest);
                    R_best = R_vals(ixBest);
                    B_best = B_vals(iyBest);

                    % Plot
                    h = surf(RR, BB, Z);
                    set(h, 'EdgeColor','none');
                    shading flat; hold on
                    plot(R_best, B_best, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
                    title(sprintf('S1=%d, S2=%d', S1_vals(ii), S2_vals(jj)), 'FontWeight','normal');
                    box on; colorbar

                end
            end
        end

        %% -------------------- BEST-OVERALL (no bootstrap weighting) --------------------
        switch upper(metric)
            case 'LL'
                scoreAll = sum(Lmat, 1);               % 1 x Ncols (bigger = better)
                [~, bestColGlobal] = max(scoreAll);
            case 'SSD'
                scoreAll = sum(Diff2, 1);              % 1 x Ncols (smaller = better)
                [~, bestColGlobal] = min(scoreAll);
        end

        pkBestGlobal = Pk(:, bestColGlobal);          % Y x 1 (full aligned years)
        bestR  = R_all(bestColGlobal);
        bestB  = B_all(bestColGlobal);
        bestS1 = S1_all(bestColGlobal);
        bestS2 = S2_all(bestColGlobal);
        eval(['save CreatedData/BestBootstrapFits_' metric '_IncludeRevisions_' num2str(IncludeRevisions)])
        fprintf('Best overall (no bootstrap): R=%.4f, B=%.4f, Span=[%g %g]\n', ...
            bestR, bestB, bestS1, bestS2);

        %% -------------------- BOOTSTRAP: BEST COLUMN PER RESAMPLE --------------------
        Pred  = nan(B, Y);                      % store best-fit lines (B x Y)
        bestIdx = zeros(B,1);                   % winning column per bootstrap

        for bs = 1:B
            % sample Yf rows with replacement from the fitting rows
            w = accumarray(randi(Yf, [Yf,1]), 1, [Yf,1]);   % multiplicities (Yf x 1)

            switch upper(metric)
                case 'LL'
                    % weighted LL per column
                    scores = (w.' * Lmat);                  % 1 x Ncols
                    [~, bestCol] = max(scores);
                case 'SSD'
                    % weighted SSD per column
                    scores = (w.' * Diff2);                 % 1 x Ncols
                    [~, bestCol] = min(scores);
            end

            bestIdx(bs) = bestCol;
            Pred(bs,:)  = Pk(:, bestCol).';                 % 1 x Y
        end

        %% -------------------- POSTERIOR HEAT (PATCHES) + QUANTILE ENVELOPE --------------------

        figure(2); clf; 
        subplot('position',[0.07 0.09 0.93 0.88]); hold on
        DIM = [0 0 35 12].*0.9;
        FS = 16;
        set(gca,'fontsize',FS-2); box on
        
        % Colour scheme
        Col = slanCM('ice',20);
        NetworkCL = Col([17 14 9 5],:);
        LineCol = Col(17,:);
        MarkCol = Col(9,:);
        gridCol = Col(18,:);
        FontCol = Col(7,:);

        % Posterior parameters for each aligned year
        aPost = s + 1;
        bPost = (n - s) + 1;

        % Y-axis grid (100 bins)
        yEdges   = linspace(0,1,101);
        Ny       = numel(yEdges)-1;

        % Rectangle width on the year axis (90% of min step)
        uYears = unique(years);
        dx = 0.8;
        if numel(uYears) >= 2
            dx = 0.9 * min(diff(uYears));
        end

        % Posterior bin mass for each year
        T  = numel(years);
        mass = zeros(Ny, T);
        for k = 1:T
            mass(:,k) = betacdf(yEdges(2:end), aPost(k), bPost(k)) - ...
                betacdf(yEdges(1:end-1), aPost(k), bPost(k));
        end

        plot_years = 3:length(years)-2;

        for p = 1:length(plot_years)
            this_n = n_u(plot_years(p));
            this_s = s_u(plot_years(p));
            this_mle = this_s./this_n;
            this_y = years_obs(plot_years(p));

            [pHat, CI] = binofit(this_s, this_n, 1-0.05); % CI is Clopper–Pearson

            pp = patch(this_y+[0 1 0 -1 0].*0.3, ...
                [CI(1) this_mle CI(2) this_mle CI(1)], MarkCol);
            set(pp,'facealpha',0.33,'edgecolor','none')

            pre(p) = pkBestGlobal(plot_years(p));
            obs(p,:) = [this_mle CI];
        end

        % % Best-overall single fit (blue)
        % plot(years(plot_years), pkBestGlobal(plot_years),...
        %     '-','color',FontCol, 'LineWidth', 2.25);
        % plot(years_obs(1:end-1), p_u(1:end-1), '.','markersize',11,'color',MarkCol)

        % ----- Best-overall single fit (blue, gently curved) -----
        x_model = years(plot_years);
        y_model = pkBestGlobal(plot_years);
        xx = linspace(min(x_model), max(x_model), 400);
        yy  = interp1(x_model, y_model, xx, 'pchip');
        plot(xx, yy, '-', 'Color', FontCol, 'LineWidth', 2.25);
        plot(x_model, y_model, 'd', 'Color', FontCol, 'markersize', 4, 'MarkerFaceColor', FontCol);

        % Observations
        plot(years_obs(1:end-1), p_u(1:end-1), '.', 'MarkerSize', 11, 'Color', MarkCol);

        TY = 2008; Th = 0.15;
        plot(TY+[0 4],[1 1].*Th,'-', 'Color', FontCol, 'LineWidth', 3);
        text(TY+4.5,Th,'FJ Model prediction','fontsize',FS-4, 'Color', FontCol);

        Th = 0.075;
        pp = patch(TY+[0 2 4 2],Th+(0.02*[0 1 0 -1]), MarkCol);
        set(pp,'facealpha',0.33,'edgecolor','none')
        text(TY+4.5,Th,'ATCM estimate','fontsize',FS-4, 'Color', FontCol);
        plot(TY+2,Th, '.','markersize',11,'color',MarkCol)


        % Axes cosmetics
        ylim([0 1]);
        xlim([min(years)-1, max(years)+1]);
        xticks(min(years):5:max(years));
        yticks(0:0.2:1);
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
            text(xB(1)+0.75,yy,sprintf(['%.0f'], yy*100),'horizontalalignment','right','verticalalignment','middle',...
                'fontsize',FS-4,'color',FontCol)
        end
        text(xB(1)-2,0.5,'Probability a Measure enters into force','rotation',90',...
            'horizontalalignment','center',...
            'fontsize',FS-2,'color',FontCol)
        axis off


        % Quick report of span distribution among winners
        if ~isempty(bestIdx)
            S1_win = S1_all(bestIdx); S2_win = S2_all(bestIdx);
            fprintf('Winners (bootstrap): median S1 = %.1f, median S2 = %.1f\n', median(S1_win), median(S2_win));
        end

    end
end