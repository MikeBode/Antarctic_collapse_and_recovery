function CreateModelTimeseries(R_grid, B_grid, NumMCRuns, seed, spanGrid)
% CREATEMODELTIMESERIES  Simulate the consensus model over the parameter grid.
%
% For each span and each combination of consensus radius (R) and log
% background influence (b), the opinion dynamics model is run NumMCRuns
% times on every year's co-authorship network. The proportion of runs that
% reach consensus is the model's predicted probability that a Measure
% enters into force that year.
%
% Inputs (all required)
%   R_grid     consensus radii
%   B_grid     log background influence values
%   NumMCRuns  Monte Carlo runs per parameter set
%   seed       random number seed
%   spanGrid   N x 2 list of [S1 S2] spans, one per row
%
% Outputs: CreatedData/Timeseries/
%   pk_cache_span_S1_S2.mat  results for one span
%   pk_flat_store.mat        all spans in one matrix, Mflat, where each
%                            column is [R; b; S1; S2; pk for each year]

if nargin < 5 || any(cellfun(@isempty, {R_grid, B_grid, NumMCRuns, seed, spanGrid}))
    error('CreateModelTimeseries:MissingArguments', ...
        'All five arguments are required: R_grid, B_grid, NumMCRuns, seed, spanGrid.');
end

inDir  = fullfile('CreatedData','AveragedNetworks');
outDir = fullfile('CreatedData','Timeseries');
if ~isfolder(outDir), mkdir(outDir); end

% Check the span list
validateattributes(spanGrid, {'numeric'}, ...
    {'2d','ncols',2,'integer','nonnegative','nonempty'}, mfilename, 'spanGrid');
if size(unique(spanGrid,'rows'),1) ~= size(spanGrid,1)
    error('CreateModelTimeseries:DuplicateSpan', 'spanGrid contains duplicate rows.');
end
nSpans = size(spanGrid, 1);

% Before any simulation, check every network file exists and that all spans
% cover the same years. Otherwise the run would fail at the very end.
yrs = cell(nSpans, 1);
for i = 1:nSpans
    fn = fullfile(inDir, sprintf('Collaboration_by_year_span_%d_%d.mat', spanGrid(i,1), spanGrid(i,2)));
    if ~isfile(fn)
        error('CreateModelTimeseries:MissingNetwork', ...
            'Network file not found: %s. Run ExtractNetworkData for this span first.', fn);
    end
    C = load(fn, 'WP_Collaboration');
    yrs{i} = cell2mat(C.WP_Collaboration(:,1)).';
end
if ~all(cellfun(@(y) isequal(y, yrs{1}), yrs))
    error('CreateModelTimeseries:YearMismatch', ...
        'Spans cover different years, so their fits would not be comparable.');
end

% Remove results from any earlier run
old = dir(fullfile(outDir, 'pk_*'));
for f = 1:numel(old)
    delete(fullfile(outDir, old(f).name));
end

% Run the model for each span
tStart = tic;
for i = 1:nSpans
    Span = spanGrid(i,:);

    % Reseed for each span so a span's results don't depend on where it
    % sits in the list.
    rng(seed);

    PrintRunConfig(R_grid, B_grid, NumMCRuns, seed, Span, ...
        sprintf('MONTE CARLO, SPAN [%d %d] (%d of %d)', Span(1), Span(2), i, nSpans));
    BuildSpanCache(Span, R_grid, B_grid, NumMCRuns, outDir, inDir);

    % Progress and estimated finish time
    elapsed = toc(tStart);
    remaining = (nSpans - i) * elapsed / i;
    fprintf('  Done %d of %d spans, elapsed %.0f s, estimated finish %s\n', ...
        i, nSpans, elapsed, char(datetime('now') + seconds(remaining)));
end

ConsolidateSpans(outDir, seed, NumMCRuns, spanGrid);
end


% ------------------------------------------------------------------------
function BuildSpanCache(Span, R_grid, B_grid, NumMCRuns, outDir, inDir)
% Runs the model over the R/b grid for one span and saves the results.

netFile = fullfile(inDir, sprintf('Collaboration_by_year_span_%d_%d', Span(1), Span(2)));
C = load(netFile, 'WP_Collaboration');
years   = cell2mat(C.WP_Collaboration(:,1)).';
AdjList = C.WP_Collaboration(:,4);
nYears  = numel(years);
N_year  = cellfun(@(A) size(A,1), AdjList).';

% Set up the output file and write into it one R value at a time
cachePath = fullfile(outDir, sprintf('pk_cache_span_%d_%d.mat', Span(1), Span(2)));
NR = numel(R_grid);
NB = numel(B_grid);
save(cachePath, '-v7.3', 'Span', 'years', 'R_grid', 'B_grid', 'NumMCRuns');
M = matfile(cachePath, 'Writable', true);
M.pk = zeros(nYears, NR, NB, 'single');

for rIdx = 1:NR
    R = R_grid(rIdx);
    pk_slice = zeros(nYears, NB, 'single');

    for bIdx = 1:NB
        BG = [exp(B_grid(bIdx)), 0];
        Consensus = false(NumMCRuns, nYears);

        for yIdx = 1:nYears
            A = AdjList{yIdx};
            if isempty(A), continue; end
            for m = 1:NumMCRuns
                Consensus(m, yIdx) = TestConsensus(N_year(yIdx), A, R, BG);
            end
        end

        pk_slice(:, bIdx) = single(mean(Consensus, 1)).';
    end

    M.pk(:, rIdx, :) = reshape(pk_slice, [nYears, 1, NB]);
end
end


% ------------------------------------------------------------------------
function ConsolidateSpans(outDir, seed, NumMCRuns, spanGrid)
% Combines the per-span files into one matrix. Files are read in spanGrid
% order so the column order is always the same.

nSpans = size(spanGrid, 1);
files = cell(nSpans, 1);
for s = 1:nSpans
    files{s} = fullfile(outDir, sprintf('pk_cache_span_%d_%d.mat', spanGrid(s,1), spanGrid(s,2)));
    assert(isfile(files{s}), 'Expected file missing: %s', files{s});
end

M0 = matfile(files{1});
years  = M0.years;
R_grid = M0.R_grid;
B_grid = M0.B_grid;
Y  = numel(years);
NR = numel(R_grid);
NB = numel(B_grid);

Mflat = zeros(4 + Y, nSpans*NR*NB, 'single');
col = 0;
for s = 1:nSpans
    Ms = matfile(files{s});
    assert(isequal(Ms.years, years), 'Years differ in %s', files{s});
    assert(isequal(Ms.R_grid, R_grid) && isequal(Ms.B_grid, B_grid), 'Grid differs in %s', files{s});

    SP    = Ms.Span;
    block = single(Ms.pk);   % Y x NR x NB

    for rIdx = 1:NR
        pk_R = reshape(block(:, rIdx, :), Y, NB);
        for bIdx = 1:NB
            col = col + 1;
            Mflat(1:4, col)   = single([R_grid(rIdx); B_grid(bIdx); SP(1); SP(2)]);
            Mflat(5:end, col) = pk_R(:, bIdx);
        end
    end
end

save(fullfile(outDir, 'pk_flat_store.mat'), 'Mflat', 'years', 'R_grid', 'B_grid', ...
     'spanGrid', 'seed', 'NumMCRuns', '-v7');
fprintf('  Saved CreatedData/Timeseries/pk_flat_store.mat (%d x %d)\n', size(Mflat,1), size(Mflat,2));
end


% ------------------------------------------------------------------------
function isConsensus = TestConsensus(N, Network, R, BG)
% One run of the opinion dynamics model. Each Party starts at a random
% point in a 2D opinion space. Over a few rounds each moves to the weighted
% average of the others' positions, weighted by co-authorship plus a small
% uniform background influence. Consensus is reached if every Party ends
% up within 2R of every other.

InitialPosition = randn(2, N);
FinalPosition   = InitialPosition;

W = Network + BG(1);

Rounds = 5;
for r = 2:Rounds
    for n = 1:N
        denom = sum(W(:,n));
        if denom == 0, denom = eps; end
        FinalPosition(:,n) = (InitialPosition * W(:,n)) / denom;
    end
    InitialPosition = FinalPosition;
end

isConsensus = max(pdist(FinalPosition')) <= 2*R;
end
