% function buildLikelihoodFromMeasures()
clear all
%
% OUTPUTS
%   Likelihood{y,1} = year (ATCM_Year)
%   Likelihood{y,2} = prob grid (row vector)
%   Likelihood{y,3} = pdf   grid (row vector)
%   Likelihood{y,4} = alpha = (# ratified in year) + 1
%   Likelihood{y,5} = beta  = (# proposed but not ratified) + 1
%
%   MeasuresInfo is a 3 x Y matrix:
%     MeasuresInfo(1,y) = ATCM_Year
%     MeasuresInfo(2,y) = total Measures proposed that year
%     MeasuresInfo(3,y) = number ratified (non-empty Adoption_Year)
%
% Notes:
% - Only rows with Type=="Measure" (case-insensitive) are used.
% - Empty/NaN in Adoption_Year counts as NOT ratified.

for DropSATCM = 0:1

    T = readtable('RawData/MeasureCorpus.csv');

    % ---------- Normalise/clean ----------
    % Ensure Type is string (lowercased for filtering)
    TypeStr = lower(strtrim(string(T.Type)));
    AdoptYear = T.Adoption_Year;
    ATCMYear  = T.ATCM_Year;
    MeetingTypeStr = T.Meeting_Type;

    % ---------- Filter to Measures and Recommendations only ----------
    isMeasure = contains(TypeStr, 'measure', 'IgnoreCase', true) | ...
        contains(TypeStr, 'recommendation', 'IgnoreCase', true);
    T = T(isMeasure, :);
    AdoptYear = AdoptYear(isMeasure);
    ATCMYear  = ATCMYear(isMeasure);

    if DropSATCM == 1 % ---------- Drop SATCMs ----------
        is_SATCM = find(strcmp(T.Meeting_Type,'ATCM'));
        T = T(is_SATCM, :);
        AdoptYear = AdoptYear(is_SATCM);
        ATCMYear  = ATCMYear(is_SATCM);
    end

    % Drop rows without an ATCM year
    good = ~isnan(ATCMYear);
    T = T(good, :);
    AdoptYear = AdoptYear(good);
    ATCMYear  = ATCMYear(good);

    % ---------- Group by ATCM_Year ----------
    years = unique(ATCMYear(:));
    Y = numel(years);

    Totals   = zeros(1, Y);
    Ratified = zeros(1, Y);

    for i = 1:Y
        y = years(i);
        idx = (ATCMYear == y);
        Totals(i)   = sum(idx);                 % total measures proposed that year
        Ratified(i) = sum(idx & ~isnan(AdoptYear));  % ratified = Adoption_Year present
    end

    % ---------- Build MeasuresInfo ----------
    MeasuresInfo = zeros(3, Y);
    MeasuresInfo(1,:) = years;           % year
    MeasuresInfo(2,:) = Totals;          % total proposed
    MeasuresInfo(3,:) = Ratified;        % ratified

    % ---------- Build Likelihood cell ----------
    % Posterior for p(year) with Beta(1+success, 1+failure)
    Likelihood = cell(Y, 5);

    % A smooth probability grid is handy for plotting
    prob = linspace(0, 1, 1001);  % same for all rows

    for i = 1:Y
        succ = Ratified(i);
        fail = Totals(i) - Ratified(i);
        alpha = succ + 1;
        beta  = fail + 1;

        pdfVals = betapdf(prob, alpha, beta);

        Likelihood{i,1} = years(i);
        Likelihood{i,2} = prob;
        Likelihood{i,3} = pdfVals;
        Likelihood{i,4} = alpha;
        Likelihood{i,5} = beta;
    end

    % ---------- Save to disk ----------
    if DropSATCM == 1
        save CreatedData/MeasureData_NoSATCM Likelihood MeasuresInfo
    elseif DropSATCM == 0
        save CreatedData/MeasureData Likelihood MeasuresInfo
    end

    %% ================= Violin plot of per-year posteriors =================
    % Requires: Likelihood{y,1}=year; {y,4}=alpha; {y,5}=beta

    % ---- Extract and order years / Beta params
    years  = cell2mat(Likelihood(:,1));
    alphaV = cell2mat(Likelihood(:,4));
    betaV  = cell2mat(Likelihood(:,5));

    % Keep unique years in ascending order (handles any duplicates/unsorted rows)
    [yearsU, ord] = sort(years(:));
    alphaV = alphaV(ord);
    betaV  = betaV(ord);

    % ---- Compute Beta pdfs on a common grid
    pGrid = linspace(0,1,400);        % probability axis (vertical)
    nY    = numel(yearsU);
    pdfMat = zeros(numel(pGrid), nY);

    for i = 1:nY
        pdfMat(:,i) = betapdf(pGrid, alphaV(i), betaV(i));
    end

    % ---- Scale widths so the widest violin is width = maxWidth
    widthMat = pdfMat./repmat(max(pdfMat),size(pdfMat,1),1)./2;

    % ---- Plot
    figure(2); clf, hold on
    CL = [0.2 0.5 0.9];
    for i = 1:nY-1
        xCenter = i;

        % Build a symmetric polygon around xCenter
        leftX  = xCenter - widthMat(:,i).';
        rightX = xCenter + widthMat(:,i).';
        yVals  = pGrid;

        % Violin body
        patch([leftX, fliplr(rightX)], [yVals, fliplr(yVals)], ...
            CL, 'FaceAlpha', 0.35, 'EdgeColor', 'none');

        % Optional: a faint center line for reference
        plot([xCenter xCenter], [0 1], ':', 'Color', [0 0 0 0.15]);

    end

    % Cosmetics
    xlim([0.5, nY+0.5]);
    ylim([0, 1]);
    xticks(1:nY);
    xticklabels(string(yearsU));
    xtickangle(45);
    ylabel('Probability a measure is ratified');
    xlabel('ATCM year');
    box on;

    % if DropSATCM == 1
    %     Make_TIFF('Figures/Figure1_ObservedMeasures_NoSATCM.tiff',[0 0 30 12])
    % elseif DropSATCM == 0
    %     Make_TIFF('Figures/Figure1_ObservedMeasures.tiff',[0 0 30 12])
    % end
end
