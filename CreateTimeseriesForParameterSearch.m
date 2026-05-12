function CreateTimeseriesForParameters()
% CREATE TIMESERIES FOR PARAMETERS
% Precompute pk time series over a 4D grid:
% Span = [S1 S2] (window in ATCM counts; endpoints excluded)
% R_grid
% B_grid -> uniform off-diagonal background = exp(B1)
%
% For each Span, writes one cache file:
% Data/Timeseries_NoRevisions/pk_cache_span_S1_S2.mat
% containing:
% - years : 1 x Y
% - Span : [S1 S2]
% - R_grid : 1 x NR
% - B1_grid : 1 x NB1
% - B2_grid : 1 x NB2
% - NumMCRuns : scalar
% - pk : Y x NR x NB1 x NB2 (single)

NoSpan = 0;
for IncludeRevisions = 1 %0:1

    % Set the load and save directories
    if IncludeRevisions == 0
        inDir = 'CreatedData/AveragedNetworks_NoRevisions';
        outDir = 'CreatedData/Timeseries_NoRevisions';
    elseif IncludeRevisions == 1
        inDir = 'CreatedData/AveragedNetworks_Revisions';
        if NoSpan == 0
            outDir = 'CreatedData/Timeseries_Revisions';
        elseif NoSpan == 1
            outDir = 'CreatedData/Timeseries_Revisions/NoSpan';
        end
    end

    % ------------ Parameters ------------
    NumMCRuns = 500; % Monte Carlo per parameter set
    GRD = 50;

    if NoSpan == 1
        R_grid = linspace(0.35, 0.65, GRD); % consensus radius
        B_grid = linspace(-6, -3, GRD); % log background (off-diagonal)
    elseif NoSpan == 0
        R_grid = linspace(0.125, 0.3, GRD); % consensus radius
        B_grid = linspace(-12, -8, GRD); % log background (off-diagonal)
    end

    % ------------ Clean up old cache files ------------
    oldFiles = dir(fullfile(outDir, 'pk_*'));
    if ~isempty(oldFiles)
        fprintf('Deleting %d existing cache files in %s ...\n', numel(oldFiles), outDir);
        for f = 1:numel(oldFiles)
            delete(fullfile(outDir, oldFiles(f).name));
        end
    end

    % ------------ Span values and progress tracking ------------
    if NoSpan == 1
        S1_vals = 0;
        S2_vals = 0;
    elseif NoSpan == 0
        S1_vals = 3;
        S2_vals = 1;
    end
    totalSpans = numel(S1_vals) * numel(S2_vals);
    spanCounter = 0; tStartAll = tic;
    fprintf('Beginning timeseries generation for %d spans (%s)\n\n',totalSpans, datestr(datetime('now')));

    % ------------ Main loop ------------
    for s1 = S1_vals
        for s2 = S2_vals
            Span = [s1 s2];
            fprintf('=== Building pk cache for Span = [%d %d] ===\n', s1, s2);
            cachePath = build_pk_cache(Span, R_grid, B_grid, NumMCRuns, outDir, inDir);

            % ---- Progress and ETA ----
            spanCounter = spanCounter + 1;
            elapsed = toc(tStartAll);
            avgPer = elapsed / spanCounter;
            remaining = (totalSpans - spanCounter) * avgPer;
            estFinish = datetime('now') + seconds(remaining);

            fprintf('Progress: %d/%d spans | elapsed %.1fs | avg/span %.1fs | ETA: %s\n\n', spanCounter, totalSpans, elapsed, avgPer, datestr(estFinish));
        end
    end

    fprintf('All spans complete. Total runtime: %.1f s (%s)\n', toc(tStartAll), datestr(datetime('now')));

    ConsolidateFits(outDir)
end
end

% ======================================================================
% Helper: Consolidate fits into single flat packed dataset
% ======================================================================
function ConsolidateFits(outDir)
% BuildFlatPkStore.m
% Consolidate all pk_cache_span_*.mat into ONE flat matrix:
%   M = [ R;
%         B;
%         S1;
%         S2;
%         pk(1);
%         pk(2);
%         ...
%         pk(Y) ]
% Each column is a parameter set (R,B,S1,S2) with its pk time series.

files = dir(fullfile(outDir, 'pk_cache_span_*.mat'));

% Read grids & years from first file (sanity-check consistency)
M0 = matfile(fullfile(files(1).folder, files(1).name));
years  = M0.years;
R_grid = M0.R_grid;
B_grid = M0.B_grid;
Y  = numel(years);
NR = numel(R_grid);
NB = numel(B_grid);

% Count total columns we will produce = sum over spans of (NR*NB)
numSpans = numel(files);
totalCols = numSpans * NR * NB;

% Preallocate final matrix: 4 param rows + Y series rows
Mflat = zeros(4 + Y, totalCols, 'single');
col = 0;

for s = 1:numSpans
    Ms = matfile(fullfile(files(s).folder, files(s).name));
    % Sanity checks
    assert(isequal(Ms.years, years),  'Years mismatch in %s', files(s).name);
    assert(isequal(Ms.R_grid, R_grid) && isequal(Ms.B_grid, B_grid), ...
        'Grid mismatch in %s', files(s).name);

    SP = Ms.Span;                 % [S1 S2]
    S1 = single(SP(1));
    S2 = single(SP(2));

    % Read block: pk is Y x NR x NB
    block = Ms.pk; % loads the whole block once (fine; single in your cache or cast)
    block = single(block);

    % Iterate R,B grid and fill columns
    for rIdx = 1:NR
        Rval = single(R_grid(rIdx));
        pk_R = block(:, rIdx, :);                % Y x 1 x NB
        pk_R = reshape(pk_R, Y, NB);             % Y x NB

        % Fill NB columns for this rIdx
        for bIdx = 1:NB
            col = col + 1;
            Bval = single(B_grid(bIdx));

            Mflat(1, col) = Rval;
            Mflat(2, col) = Bval;
            Mflat(3, col) = S1;
            Mflat(4, col) = S2;
            Mflat(5:end, col) = pk_R(:, bIdx);
        end
    end

    if mod(s,5)==0 || s==numSpans
        fprintf('Packed span %d/%d (S=[%g %g])\n', s, numSpans, S1, S2);
    end
end

% Save once
outMat = fullfile(outDir,'pk_flat_store.mat');
save(outMat, 'Mflat', 'years', 'R_grid', 'B_grid', '-v7'); % v7 is fine here
fprintf('Saved flat store: %s  (size: %dx%d)\n', outMat, size(Mflat,1), size(Mflat,2));
end

% ======================================================================
% Helper: build cache for a *single* Span (no parfor writes to matfile)
% ======================================================================
function cachePath = build_pk_cache(Span, R_grid, B_grid, NumMCRuns, outDir, inDir)

S1 = Span(1); S2 = Span(2);

% ----------------- Build/Load annual networks once for this Span -----------------
eval(['C = load(''' inDir '/Collaboration_by_year_span_' num2str(S1) '_' num2str(S2) ''', ''WP_Collaboration'');'])

years = cell2mat(C.WP_Collaboration(:,1)).';
AdjList = C.WP_Collaboration(:,4);
NumATCMs = numel(years);

% Precompute matrix sizes per year (for speed)
N_year = zeros(1, NumATCMs);
for y = 1:NumATCMs
    A = AdjList{y};
    N_year(y) = size(A,1);
end

% ----------------- Prepare on-disk cache (header + allocate) -----------------
cachePath = fullfile(outDir, sprintf('pk_cache_span_%d_%d.mat', S1, S2));
NR = numel(R_grid);
NB = numel(B_grid);

save(cachePath, '-v7.3', 'Span', 'years', 'R_grid', 'B_grid', 'NumMCRuns');
M = matfile(cachePath, 'Writable', true);
M.pk = zeros(NumATCMs, NR, NB, 'single'); % allocate Y × NR × NB

% ----------------- Compute per-R slabs (no parallel to simplify I/O) -----------------
for rIdx = 1:NR
    R = R_grid(rIdx);
    pk_slice = zeros(NumATCMs, NB, 'single');

    for bIdx = 1:NB
        BG = [exp(B_grid(bIdx)), 0];
        Consensus = false(NumMCRuns, NumATCMs);
        disp([R BG])

        % Loop over ATCM years
        for yIdx = 1:NumATCMs
            A = AdjList{yIdx};
            if isempty(A), continue; end
            N = N_year(yIdx);

            % Monte Carlo sampling over initial conditions
            for m = 1:NumMCRuns
                Consensus(m, yIdx) = TestConsensus(N, A, R, BG);
            end
        end

        % Probability per year
        pk_slice(:, bIdx) = single(mean(Consensus, 1)).';
    end

    % Store the slab (reshape for matfile)
    M.pk(:, rIdx, :) = reshape(pk_slice, [NumATCMs, 1, NB]);
end
end


% ======================================================================
% Local copy of TestConsensus (same logic as your model)
% ======================================================================
function isConsensus = TestConsensus(N, Network, R, BG)

T = 2; % 2D opinion space
InitialPosition = randn(T, N);
FinalPosition = InitialPosition;

% Construct weighted network
W = Network;
W = W + BG(1); % uniform off-diagonal influence
% W(1:N+1:end) = BG(2); % diagonal self-weight (stubbornness)

Rounds = 5;
for r = 2:Rounds
    for n = 1:N
        denom = sum(W(:,n));
        if denom == 0, denom = eps; end
        FinalPosition(:,n) = (InitialPosition * W(:,n)) / denom;
    end
    InitialPosition = FinalPosition;
end

% Measure consensus
D = pdist(FinalPosition'); % pairwise distances
isConsensus = max(D) <= 2*R; % success if all within 2R
end