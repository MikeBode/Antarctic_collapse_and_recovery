function FitModel(fitYears, nBoot)
% FITMODEL  Find the model run that best matches the observed data.
%
% Compares every column of the simulated output against the observed
% proportion of Measures entering into force each year, using the sum of
% squared differences (SSR) over the fit window. Years with no Measures
% adopted are left out. A bootstrap over years shows how stable the best
% fit is.
%
% Inputs
%   fitYears  [first last] years used in the fit
%   nBoot     number of bootstrap resamples
%
% Output: CreatedData/BestFit.mat

M = load(fullfile('CreatedData','MeasureData.mat'), 'Likelihood');
S = load(fullfile('CreatedData','Timeseries','pk_flat_store.mat'), ...
    'Mflat', 'years', 'seed', 'R_grid', 'B_grid', 'NumMCRuns');
Mflat = S.Mflat;
seed  = S.seed;
years_flat = S.years(:);

PrintRunConfig(S.R_grid, S.B_grid, S.NumMCRuns, seed, [], 'FIT CONFIGURATION');

% Parameter values and simulated time series for each column
R_all  = double(Mflat(1,:));
B_all  = double(Mflat(2,:));
S1_all = double(Mflat(3,:));
S2_all = double(Mflat(4,:));
Pk_all = double(Mflat(5:end,:));

% Observed counts and proportions per year
Likelihood = M.Likelihood;
y_raw = cell2mat(Likelihood(:,1));
a_raw = cell2mat(Likelihood(:,4));
b_raw = cell2mat(Likelihood(:,5));
[years_obs, ~, g] = unique(y_raw(:));
s_u = accumarray(g, a_raw(:) - 1, [], @sum);           % entered into force
n_u = accumarray(g, a_raw(:) + b_raw(:) - 2, [], @sum); % adopted
p_u = s_u ./ max(n_u, 1);

% Line up model and observed years
[years, iFlat, iObs] = intersect(years_flat, years_obs, 'stable');
if isempty(years)
    error('No years in common between the model output and the observations.');
end
Pk   = Pk_all(iFlat, :);
s    = double(s_u(iObs));
n    = double(n_u(iObs));
pobs = double(p_u(iObs));

% Years used in the fit
fitMask = years >= fitYears(1) & years <= fitYears(2) & n > 0;
Yf = sum(fitMask);

% Squared differences for every column, and the overall best fit
Diff2 = (Pk(fitMask,:) - pobs(fitMask)).^2;
scoreAll = sum(Diff2, 1);
[~, bestColGlobal] = min(scoreAll);

pkBestGlobal = Pk(:, bestColGlobal);
bestR  = R_all(bestColGlobal);
bestB  = B_all(bestColGlobal);
bestS1 = S1_all(bestColGlobal);
bestS2 = S2_all(bestColGlobal);
fprintf('  Best fit: R = %.4f, b = %.4f, span = [%g %g]\n', bestR, bestB, bestS1, bestS2);

% Bootstrap: resample the fit years with replacement and refit each time
B = nBoot;
rng(seed);
Pred    = nan(B, numel(years));
bestIdx = zeros(B, 1);
for bs = 1:B
    w = accumarray(randi(Yf, [Yf,1]), 1, [Yf,1]);
    [~, bestCol] = min(w.' * Diff2);
    bestIdx(bs) = bestCol;
    Pred(bs,:)  = Pk(:, bestCol).';
end

save(fullfile('CreatedData','BestFit.mat'), 'seed', ...
    'bestR', 'bestB', 'bestS1', 'bestS2', 'bestColGlobal', ...
    'pkBestGlobal', 'years', 'scoreAll', 'R_all', 'B_all', 'S1_all', 'S2_all', ...
    'fitYears', 'fitMask', 'Yf', 's', 'n', 'pobs', ...
    'years_obs', 's_u', 'n_u', 'B', 'bestIdx', 'Pred');
fprintf('  Saved CreatedData/BestFit.mat\n');
end
