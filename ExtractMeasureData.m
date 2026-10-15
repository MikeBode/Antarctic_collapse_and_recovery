function ExtractMeasureData()
% EXTRACTMEASUREDATA  Annual counts of Measures adopted and entered into force.
%
% Reads RawData/MeasureCorpus.csv, keeps Measures and Recommendations
% adopted at ATCMs, and counts for each meeting year how many were adopted
% and how many have since entered into force. A Measure with no value in
% Adoption_Year is treated as not in force.
%
% Output: CreatedData/MeasureData.mat
%   MeasuresInfo  3 x Y matrix: [year; number adopted; number in force]
%   Likelihood    Y x 5 cell: {year, probability grid, Beta pdf, alpha, beta}
%                 where alpha = in force + 1 and beta = not in force + 1

T = readtable(fullfile('RawData','MeasureCorpus.csv'));

% Keep Measures and Recommendations only
TypeStr = lower(strtrim(string(T.Type)));
keep = contains(TypeStr, 'measure') | contains(TypeStr, 'recommendation');
T = T(keep, :);

% Keep ordinary ATCMs (drops Special ATCMs)
T = T(strcmp(T.Meeting_Type, 'ATCM'), :);

AdoptYear = T.Adoption_Year;
ATCMYear  = T.ATCM_Year;

% Count adopted and in force for each meeting year
years = unique(ATCMYear(:));
Y = numel(years);
Totals   = zeros(1, Y);
Ratified = zeros(1, Y);
for i = 1:Y
    idx = (ATCMYear == years(i));
    Totals(i)   = sum(idx);
    Ratified(i) = sum(idx & ~isnan(AdoptYear));
end

MeasuresInfo = [years(:).'; Totals; Ratified];

% Beta posterior for each year's probability of entering into force
prob = linspace(0, 1, 1001);
Likelihood = cell(Y, 5);
for i = 1:Y
    alpha = Ratified(i) + 1;
    beta  = Totals(i) - Ratified(i) + 1;
    Likelihood{i,1} = years(i);
    Likelihood{i,2} = prob;
    Likelihood{i,3} = betapdf(prob, alpha, beta);
    Likelihood{i,4} = alpha;
    Likelihood{i,5} = beta;
end

save(fullfile('CreatedData','MeasureData.mat'), 'Likelihood', 'MeasuresInfo');
fprintf('  Saved CreatedData/MeasureData.mat (%d meeting years)\n', Y);
end
