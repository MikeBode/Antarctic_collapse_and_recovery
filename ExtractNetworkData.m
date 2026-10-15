function ExtractNetworkData(spanList)
% EXTRACTNETWORKDATA  Co-authorship networks between Consultative Parties.
%
% Two Parties are linked when they co-author a Working Paper. For each ATCM
% the network counts co-authored papers over a window of meetings: S1
% meetings before and S2 meetings after the focal one. Only Parties that
% were Consultative by that year are kept.
%
% Input
%   spanList  N x 2 list of [S1 S2] spans to build
%
% Output: CreatedData/AveragedNetworks/
%   CoauthorshipData.mat                  cleaned author and party lists
%   Collaboration_by_year_span_S1_S2.mat  one file per span, containing
%     WP_Collaboration           networks over the span window
%     WP_Collaboration_thisyear  networks for each meeting alone
%   Each row of these is {year, number of CPs, CP names,
%   CP-only matrix, full matrix, total weight, number of papers}.

outDir = fullfile('CreatedData','AveragedNetworks');

% Working Papers and their authors. The sheet labels type as both 'WP' and
% 'wp' (the lowercase rows are University of Tasmania archive papers), so
% the match is case-insensitive.
[D, T] = xlsread(fullfile('RawData','ATCM_IPWP.xlsx'));
T = T(2:end, :);
Author = T(:,8);
Type   = T(:,6);
Year   = D(:,1);

isWP = strcmpi(Type, 'WP');
Author_wp = Author(isWP);
Year_wp   = Year(isWP);

% Consultative Parties and the year each gained Consultative status
[Dc, Tc] = xlsread(fullfile('RawData','Country_joined.xlsx'));
Parties = Tc(:,1);
isCP = ~isnan(Dc(:,2));
Consult     = Parties(isCP);
YearConsult = Dc(isCP, 2);

save(fullfile(outDir, 'CoauthorshipData'), 'Author_wp', 'Year_wp', 'Consult', 'YearConsult');

% Build the networks for each requested span
for s = 1:size(spanList, 1)
    fprintf('  Building networks for span [%d %d]\n', spanList(s,1), spanList(s,2));
    BuildSpanNetworks(spanList(s,:), Author_wp, Year_wp, Consult, YearConsult, outDir);
end
end


% ------------------------------------------------------------------------
function BuildSpanNetworks(Span, Author_wp, Year_wp, Consult, YearConsult, outDir)
% Builds and saves the networks for one [S1 S2] span.

ATCM_years = unique(Year_wp);
[~, idx_wp] = ismember(Year_wp, ATCM_years);   % meeting number for each paper
nYears = numel(ATCM_years);

WP_Collaboration          = cell(nYears, 7);
WP_Collaboration_thisyear = cell(nYears, 7);

for y = 1:nYears

    % Parties that were Consultative by this year
    ThisCP = find(YearConsult <= ATCM_years(y));

    % Papers within the span window around this meeting
    lowIdx  = max(1, y - Span(1));
    highIdx = min(nYears, y + Span(2));
    F = find(idx_wp >= lowIdx & idx_wp <= highIdx);
    A = CountCoauthorship(Author_wp(F), Consult);

    WP_Collaboration{y,1} = ATCM_years(y);
    WP_Collaboration{y,2} = numel(ThisCP);
    WP_Collaboration{y,3} = Consult(ThisCP);
    WP_Collaboration{y,4} = A(ThisCP, ThisCP);
    WP_Collaboration{y,5} = A;
    WP_Collaboration{y,6} = sum(A(:));
    WP_Collaboration{y,7} = numel(F);

    % Papers from this meeting only
    F = find(Year_wp == ATCM_years(y));
    A = CountCoauthorship(Author_wp(F), Consult);

    WP_Collaboration_thisyear{y,1} = ATCM_years(y);
    WP_Collaboration_thisyear{y,2} = numel(ThisCP);
    WP_Collaboration_thisyear{y,3} = Consult(ThisCP);
    WP_Collaboration_thisyear{y,4} = A(ThisCP, ThisCP);
    WP_Collaboration_thisyear{y,5} = A;
    WP_Collaboration_thisyear{y,6} = sum(A(:));
    WP_Collaboration_thisyear{y,7} = numel(F);
end

outName = fullfile(outDir, sprintf('Collaboration_by_year_span_%d_%d', Span(1), Span(2)));
save(outName, 'WP_Collaboration', 'WP_Collaboration_thisyear');
end


% ------------------------------------------------------------------------
function A = CountCoauthorship(Authors, Consult)
% Symmetric matrix of co-authorship counts between Consultative Parties.
% A paper only counts if its author string lists more than one submitter
% (contains a comma). Each pair of Parties on the paper gets +1.

n = numel(Consult);
A = zeros(n);
for w = 1:numel(Authors)
    str = Authors{w};
    if ~contains(str, ','), continue; end
    idx = find(cellfun(@(c) contains(str, c), Consult));
    A(idx, idx) = A(idx, idx) + 1;
end
A(1:n+1:end) = 0;   % no self-links
end
