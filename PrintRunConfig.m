function PrintRunConfig(R_grid, B_grid, NumMCRuns, seed, Span, label)
% PRINTRUNCONFIG  Prints the settings a step is running with.
%
%   R_grid     consensus radii
%   B_grid     log background influence values
%   NumMCRuns  Monte Carlo runs per parameter set
%   seed       random number seed
%   Span       [S1 S2], or [] to leave out
%   label      heading for the printout

fprintf('\n==== %s ====\n', label);
fprintf('  seed       : %d\n', seed);
fprintf('  R grid     : %d points, %.6f to %.6f\n', numel(R_grid), min(R_grid), max(R_grid));
fprintf('  b grid     : %d points, %.6f to %.6f\n', numel(B_grid), min(B_grid), max(B_grid));
fprintf('  grid size  : %d parameter sets\n', numel(R_grid)*numel(B_grid));
fprintf('  MC runs    : %d\n', NumMCRuns);
if ~isempty(Span)
    fprintf('  span       : [%d %d]\n', Span(1), Span(2));
end
end
