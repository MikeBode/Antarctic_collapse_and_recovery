function CreateAllNetworksSpans()
% Runs CreateAnnualNetworks for Span = -1, then all integer [i j] with i,j=0..7
% Each run’s outputs are saved as: Data/Collab_by_year_span_<i>_<j>.mat

% --- 1) Prime the cache (reads and caches core inputs only)

for IncludeRevisions = 0:1

    % Set the save directory
    if IncludeRevisions == 0
        DIR = 'CreatedData/AveragedNetworks_NoRevisions';
    elseif IncludeRevisions == 1
        DIR = 'CreatedData/AveragedNetworks_Revisions';
    end

    % Initialise the dataset
    CreateAnnualNetworks(-1,IncludeRevisions,DIR);

    % --- 2) Loop over all integer spans [i j], i,j in 0..7
    for i = 0:7
        for j = 0:7
            fprintf('Running Span = [%d %d] ...\n', i, j); tspan = tic;
            % Run the network construction
            CreateAnnualNetworks([i j],IncludeRevisions,DIR);
        end
    end
end
end


function CreateAnnualNetworks(Span,IncludeRevisions,DIR)
% CREATEANNUALNETWORKS builds co-authorship adjacency matrices for ATCM
% Working Papers (WPs), both (i) within a multi-ATCM "span" around each year
% and (ii) for each ATCM alone.
%
% INPUT
%   Span : either -1 or a 1x2 vector [pastATCMs, futureATCMs]
%          - If Span == -1, this branch ONLY reads and caches the core data
%            (authors, types, years, consultative parties) to MAT files.
%            No collaboration matrices are computed here.
%          - Otherwise, it expects the cached data to exist and computes:
%              * WP_Collaboration          — matrices built from a moving window
%                of ATCMs around each focal ATCM (excludes the end indices).
%              * WP_Collaboration_thisyear — matrices built from that year only.
%
%   IncludeRevisions = 1 (yes) or 0 (no)
%
%   DIR = saved file directory
%
% -------------------------------------------------------------------------

if Span == -1

    % =========================
    % 1) LOAD & PREP THE INPUTS
    % =========================

    % Load main table of papers (authors/types/years).
    % xlsread returns numeric (D) and cell (T) parts. 'Headers' are first row.
    [D,T] = xlsread('RawData/ATCM_IPWP.xlsx');
    Headers = T(1,:);             % cell array of column names (not used later)
    T = T(2:end,:);               % drop header row from cell data

    Author   = T(:,8);              % author/affiliation string (cellstr)
    Type     = T(:,6);              % 'wp', 'ip', 'sp', etc. (cellstr)
    Year     = D(:,1);              % numeric year vector
    Meeting  = D(:,2);              % ATCM Meeting
    Revision = D(:,9);              % ATCM Meeting

    % Keep only Working Papers (wp) for later network construction
    F = strcmpi(Type,'WP');
    Author_wp   = Author(F,:);      % authors for WPs (cell array)
    Year_wp     = Year(F);          % years for WPs   (numeric)
    Meeting_wp  = Meeting(F);       % years for WPs   (numeric)
    Revision_wp = Revision(F);

    if IncludeRevisions == 0
        % Remove duplicates that are coming through from revisions
        F = find(Revision_wp == 0);
        Author_wp   = Author_wp(F,:);      % authors for WPs (cell array)
        Year_wp     = Year_wp(F);          % years for WPs   (numeric)
        Meeting_wp  = Meeting_wp(F);       % years for WPs   (numeric)
        Revision_wp = Revision_wp(F);
    end

    % Load consultative party metadata:
    %   Parties list in Tc(:,1)
    %   Dc(:,2) is the year the Party becomes Consultative (NaN = non-consultative)
    [Dc,Tc] = xlsread('RawData/Country_joined.xlsx');
    Parties     = Tc(:,1);                         % all Parties (names)
    NonConsult  = Parties(isnan(Dc(:,2)));         % never consultative
    Consult     = Parties(isnan(Dc(:,2))==0);      % consultative at some point
    YearConsult = Dc(isnan(Dc(:,2))==0,2);         % year each became consultative

    % Save the core inputs to reuse during matrix construction
    ['save ' DIR '/CoauthorshipData']
    eval(['save ' DIR '/CoauthorshipData'])

else

    % =================================================
    % 2) BUILD NETWORKS (MOVING WINDOW & YEAR-BY-YEAR)
    % =================================================

    % Load only what we need from the preprocessed cache
    % load Data/CoauthorshipData Author_wp Year_wp Consult YearConsult
    eval(['load ' DIR '/CoauthorshipData Author_wp Year_wp Consult YearConsult'])
    NumConsult = length(Consult);
    NumWPs = 0; % Running total of number of working papers
    NumCoWPs = 0; % Running total of coauthored working papers

    % Unique ATCM years present in the WP data (ascending)
    ATCM_years = unique(Year_wp);

    % Precompute the ATCM-index for each WP's year (so we can window by ATCM counts)
    [~, idx_wp] = ismember(Year_wp, ATCM_years);  % idx_wp is in 1..numel(ATCM_years)

    % Loop over each focal ATCM year
    for y = 1:length(ATCM_years)

        % -------------------------------
        % (a) Identify CPs active this year
        % -------------------------------
        % Consultative Parties are those whose consultative year <= focal year
        ThisCP          = find(YearConsult <= ATCM_years(y));
        ThisNumConsult  = length(ThisCP);
        ThisNameConsult = Consult(ThisCP);

        % -----------------------------------------------------------
        % (b) WINDOWED (BY ATCM INDEX): Select WPs strictly inside the span
        %     of ATCM meetings around the focal ATCM index y.
        %     Lower index bound:   max(1, y-Span(1))   (exclusive)
        %     Upper index bound:   min(end, y+Span(2)) (exclusive)
        %     => Endpoints excluded by '>' and '<' on ATCM indices.
        % -----------------------------------------------------------
        lowIdx  = max(1, y - Span(1));
        highIdx = min(length(ATCM_years), y + Span(2));
        F = find(idx_wp >= lowIdx & idx_wp <= highIdx);

        ThisPeriodAuthor = Author_wp(F);   % authors for WPs in the indexed span
        ThisNumWP        = length(F);      % number of WPs used for this span

        % Initialise full collaboration matrix over *all* consultative parties
        ThisCollaborationMatrix = zeros(NumConsult);

        % ----------------------------------------------------------
        % (c) Build co-authorship counts for the windowed set of WPs
        %     We parse the author string and look for CP name substrings.
        %     For each paper, we increment edges between all distinct pairs.
        % ----------------------------------------------------------
        for w = 1:length(ThisPeriodAuthor)
            NumWPs = NumWPs + 1;

            % String for this WP’s author/affiliation line
            ThisWPAuthors = ThisPeriodAuthor{w};

            % Only proceed if there’s a comma — i.e., at least 2 authors listed
            if isempty(strfind(ThisWPAuthors,',')) == 0
                NumCoWPs = NumCoWPs + 1;

                ThisAuthorNums = [];  % indices into 'Consult' that appear on this WP

                % For each consultative party name, check if it appears in the string
                for a = 1:length(Consult)
                    if strfind(ThisWPAuthors, Consult{a})
                        ThisAuthorNums = [ThisAuthorNums a];
                    end
                end

                % Make sure that there are no duplications
                ThisAuthorNums = unique(ThisAuthorNums, 'stable');

                % Increment pairwise co-authorship counts only for j > i
                for i = 1:length(ThisAuthorNums)
                    for j = i+1:length(ThisAuthorNums)
                        ThisCollaborationMatrix(ThisAuthorNums(i), ThisAuthorNums(j)) = ...
                            ThisCollaborationMatrix(ThisAuthorNums(i), ThisAuthorNums(j)) + 1;
                    end
                end
            end
        end
        % Mirror once to make the matrix symmetric (undirected)
        ThisCollaborationMatrix = ThisCollaborationMatrix + ThisCollaborationMatrix.';

        % Restrict to CPs that are consultative *by this year*
        SmallCollaborationMatrix = ThisCollaborationMatrix(ThisCP, ThisCP);

        % Store windowed results for this focal year
        WP_Collaboration{y,1} = ATCM_years(y);               % focal year
        WP_Collaboration{y,2} = ThisNumConsult;               % # consultative parties this year
        WP_Collaboration{y,3} = ThisNameConsult;              % names of CPs this year
        WP_Collaboration{y,4} = SmallCollaborationMatrix;     % windowed matrix, restricted to ThisCP
        WP_Collaboration{y,5} = ThisCollaborationMatrix;      % windowed matrix, full CP index set
        WP_Collaboration{y,6} = sum(ThisCollaborationMatrix(:)); % total co-authorship weight
        WP_Collaboration{y,7} = ThisNumWP;                    % # WPs used in window

        clearvars -except WP_Collaboration* ATCM_years y DIR ThisCP NumWPs NumCoWPs ThisNumConsult ThisNameConsult Author_wp Year_wp Consult YearConsult NumConsult Span idx_wp

        % -------------------------------------------------
        % (d) THIS-YEAR-ONLY: Select WPs exactly in year y
        % -------------------------------------------------
        F = find(Year_wp == ATCM_years(y));
        ThisNumWP        = length(F);
        ThisPeriodAuthor = Author_wp(F);

        % Build the per-year co-authorship matrix (same logic as above)
        ThisCollaborationMatrix = zeros(NumConsult);

        for w = 1:length(ThisPeriodAuthor)
            ThisWPAuthors = ThisPeriodAuthor{w};

            if isempty(strfind(ThisWPAuthors,',')) == 0
                ThisAuthorNums = [];

                for a = 1:length(Consult)
                    if strfind(ThisWPAuthors, Consult{a})
                        ThisAuthorNums = [ThisAuthorNums a]; %#ok<AGROW>
                    end
                end

                for a = 1:length(ThisAuthorNums)
                    for b = 1:length(ThisAuthorNums)
                        if a ~= b
                            ThisCollaborationMatrix(ThisAuthorNums(a), ThisAuthorNums(b)) = ...
                                ThisCollaborationMatrix(ThisAuthorNums(a), ThisAuthorNums(b)) + 1;
                        end
                    end
                end
            end
        end

        % Restrict to CPs that are consultative *by this year*
        SmallCollaborationMatrix = ThisCollaborationMatrix(ThisCP, ThisCP);

        % Store per-year results
        WP_Collaboration_thisyear{y,1} = ATCM_years(y);
        WP_Collaboration_thisyear{y,2} = ThisNumConsult;
        WP_Collaboration_thisyear{y,3} = ThisNameConsult;
        WP_Collaboration_thisyear{y,4} = SmallCollaborationMatrix;       % year-only, restricted
        WP_Collaboration_thisyear{y,5} = ThisCollaborationMatrix;        % year-only, full
        WP_Collaboration_thisyear{y,6} = sum(ThisCollaborationMatrix(:)); % total weight
        WP_Collaboration_thisyear{y,7} = ThisNumWP;                       % # WPs this year

        clearvars -except WP_Collaboration* ATCM_years DIR y NumWPs NumCoWPs ThisCP ThisNumConsult ThisNameConsult Author_wp Year_wp Consult YearConsult NumConsult Span idx_wp
    end

    eval(['save ' DIR '/Collaboration_by_year_span_' num2str(Span(1)) '_' num2str(Span(2)) ' WP_Collaboration*'])

end
end