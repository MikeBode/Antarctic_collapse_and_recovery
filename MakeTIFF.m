function MakeTIFF(FileName, Dimensions, Resolution)
% MAKETIFF  Saves the current figure as a TIFF.
%
%   FileName    output path, without extension
%   Dimensions  [left bottom width height] in centimetres
%   Resolution  optional, e.g. '-r300' (default '-r200')

if nargin < 3
    Resolution = '-r200';
end

set(gcf, 'paperunits', 'centimeters', 'paperposition', Dimensions)
set(gcf, 'renderer', 'painters')
print('-dtiff', Resolution, FileName)
end
