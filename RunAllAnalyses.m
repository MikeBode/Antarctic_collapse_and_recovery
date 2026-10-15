%% RunAllAnalyses.m
% Reproduces the analysis and figures for the paper from the raw data.
% Run this from the repository folder. It creates the output folders,
% then works through five steps:
%   1. Count Measures adopted and entered into force at each ATCM
%   2. Build the Working Paper co-authorship networks
%   3. Simulate the consensus model over the parameter grid (slow)
%   4. Fit the model to the observed data
%   5. Draw Figures 1 and 2

clear; clc;

% Run mode
%   'full'    runs every step, including the Monte Carlo (about 9 hours)
%   'figures' skips step 3 and reuses CreatedData/Timeseries/pk_flat_store.mat
RunMode = 'full';

% Model settings. These define the parameter grid used in the paper.
Settings.seed = 1;        % random number seed
Settings.R0   = 0.16;     % centre of the consensus radius (R) grid
Settings.b0   = -9.5;     % centre of the log background influence (b) grid
Settings.frac = 0.2;      % grid extends this fraction either side of the centre
Settings.n    = 25;       % grid points along each axis
Settings.D    = 2500;     % Monte Carlo runs per parameter set

[s1, s2] = ndgrid(0:3, 0:3);
Settings.spanGrid = [s1(:) s2(:)];   % spans searched by the model
Settings.spanFig1 = [0 0];           % span drawn in Figure 1

Settings.fitYears = [1964 2020];     % years used when fitting the model
Settings.nBoot    = 200;             % bootstrap resamples

R_grid = linspace(Settings.R0*(1 - Settings.frac), Settings.R0*(1 + Settings.frac), Settings.n);
b_grid = linspace(Settings.b0 - abs(Settings.b0)*Settings.frac, ...
                  Settings.b0 + abs(Settings.b0)*Settings.frac, Settings.n);

% Check the run mode and that we are in the repository folder
if ~any(strcmpi(RunMode, {'full','figures'}))
    error('RunMode must be ''full'' or ''figures''.');
end
if ~isfile(fullfile('RawData','MeasureCorpus.csv'))
    error('RawData/MeasureCorpus.csv not found. Run this from the repository folder.');
end

% Create the output folders if they are not already there
outFolders = {'CreatedData', ...
              fullfile('CreatedData','AveragedNetworks'), ...
              fullfile('CreatedData','Timeseries'), ...
              'Figures'};
for k = 1:numel(outFolders)
    if ~isfolder(outFolders{k}), mkdir(outFolders{k}); end
end

% Step 1: Measure counts
fprintf('\nStep 1 of 5: Measure data\n');
ExtractMeasureData();

% Step 2: Co-authorship networks for every span needed later
fprintf('\nStep 2 of 5: Co-authorship networks\n');
netSpans = unique([Settings.spanFig1; Settings.spanGrid], 'rows', 'stable');
ExtractNetworkData(netSpans);

% Step 3: Monte Carlo simulation, or check the stored run matches the settings
if strcmpi(RunMode, 'full')
    fprintf('\nStep 3 of 5: Monte Carlo simulation\n');
    CreateModelTimeseries(R_grid, b_grid, Settings.D, Settings.seed, Settings.spanGrid);
else
    fprintf('\nStep 3 of 5: Skipped, reusing stored model output\n');
    CheckStoredModelRun(R_grid, b_grid, Settings);
end

% Step 4: Fit the model
fprintf('\nStep 4 of 5: Model fit\n');
FitModel(Settings.fitYears, Settings.nBoot);

% Step 5: Figures
fprintf('\nStep 5 of 5: Figures\n');
Figure1_MeasuresAndNetworks(Settings.spanFig1);
Figure2_ModelledObserved();

% Short summary of the best fit, including whether it sits on the grid edge
F = load(fullfile('CreatedData','BestFit.mat'));
[~, iR] = min(abs(R_grid - F.bestR));
[~, ib] = min(abs(b_grid - F.bestB));
onEdge = iR == 1 || iR == numel(R_grid) || ib == 1 || ib == numel(b_grid);

fprintf('\nBest fit\n');
fprintf('  R    : %.6f  (grid point %d of %d)\n', F.bestR, iR, numel(R_grid));
fprintf('  b    : %.6f  (grid point %d of %d)\n', F.bestB, ib, numel(b_grid));
fprintf('  span : [%g %g]\n', F.bestS1, F.bestS2);
fprintf('  SSR  : %.6g\n', min(F.scoreAll));
if onEdge
    fprintf('  Warning: the best fit is on the edge of the R/b grid.\n');
end
fprintf('\nDone. Figures are in Figures/, data products in CreatedData/.\n');


% ------------------------------------------------------------------------
function CheckStoredModelRun(R_grid, b_grid, Settings)
% Makes sure the stored Monte Carlo output was produced with the same grid,
% seed, number of runs and spans as the settings above.

store = fullfile('CreatedData','Timeseries','pk_flat_store.mat');
if ~isfile(store)
    error('No pk_flat_store.mat found. Set RunMode = ''full'' for the first run.');
end
S = load(store, 'R_grid', 'B_grid', 'seed', 'NumMCRuns', 'spanGrid');

same = @(a,b) numel(a) == numel(b) && max(abs(a(:) - b(:))) < 1e-12;
problems = {};
if ~same(S.R_grid, R_grid),  problems{end+1} = 'R grid';          end
if ~same(S.B_grid, b_grid),  problems{end+1} = 'b grid';          end
if S.seed ~= Settings.seed,  problems{end+1} = 'seed';            end
if S.NumMCRuns ~= Settings.D, problems{end+1} = 'Monte Carlo runs'; end
if ~isequal(sortrows(double(S.spanGrid)), sortrows(double(Settings.spanGrid)))
    problems{end+1} = 'span grid';
end

if ~isempty(problems)
    error(['The stored model output does not match the current settings (%s). ' ...
           'Change the settings back or set RunMode = ''full''.'], strjoin(problems, ', '));
end
fprintf('  Stored model output matches the current settings.\n');
end
