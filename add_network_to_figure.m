function Qvs = add_network_to_figure(ThisNetwork, ThisYear, r, YS, Title, TShift, FS, LW, MKT, MKS, CL, W)
% ADD_NETWORK_TO_FIGURE
% -------------------------------------------------------------------------
% PURPOSE
%   Draw a small, circular “network glyph” for a given ATCM year at a
%   specified (x,y) location in the current axes. Nodes are placed evenly
%   on a circle; edges are drawn between node pairs with weight > 0.
%   Heavier edges get warmer colours, while light edges are drawn in grey.
%   The glyph can be rendered on a white or dark background “badge” to
%   visually separate it from the figure.
%
% INPUTS
%   ThisNetwork : [N x N double] symmetric matrix of (undirected) weights
%                 for the consultative parties in that year (zero = no edge).
%   ThisYear    : [scalar double] x-position of the glyph (often the year).
%   r           : [1 x 2 double] radii scaling [rx ry] for the node circle
%                 in x and y (controls glyph size/aspect).
%   YS          : [scalar double] y-position of the glyph (vertical center).
%   Title       : [char/string] text label drawn above the glyph
%                 (e.g., 'ATCM XLIII' or year text). Empty '' disables text.
%   TShift      : [scalar double] vertical offset for the text label (units
%                 in figure’s y axis; positive = above the glyph).
%   FS          : [scalar double] font size for the label.
%   LW          : [1 x 2 double] line widths [LW_light LW_heavy] for
%                 (1) grey/light edges and (2) coloured/heavier edges.
%   MKT         : [char] node marker type (e.g., 'k.' for small black dots).
%   MKS         : [scalar double] marker size for node markers.
%   W           : [integer] background mode:
%                   0 = no badge (transparent)
%                   1 = white circular badge behind the glyph
%                   2 = dark circular badge behind the glyph + white nodes
%                   3 = “nodes-only” on dark badge (no edges, larger white nodes)
%
% RENDERING NOTES
%   - Nodes are placed at equally spaced angles around the circle.
%   - Edges are drawn in two passes:
%       1) Grey pass for all edges with small weight (for context).
%       2) Coloured pass for edges above weight quantiles.
%   - Edge colour bins are defined by empirical quantiles of nonzero weights.
%   - If W==3, edges are suppressed and only prominent white nodes are drawn.
%
% EXAMPLE
%   % Minimal example within a figure:
%   % (Assume A is your N x N adjacency for a specific year)
%   figure; hold on; axis equal off
%   add_network_to_figure(A, 2005, [1.3 1.3], 0.75, 'ATCM XXVIII', 0.08, 10, [0.5 1.25], 'k.', 6, 1);
%
%   % Matches earlier calls you used, e.g.:
%   % add_network_to_figure(WP_Collaboration{i,4}, x(i), r, y(i), t_dates{i,1}, TShift, FS, LW, MKT, MKS, 1)
%
% DEPENDENCIES
%   - Expects a colormap function VIRIDIS on path (or replace with colormap call).
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 1) Colour binning based on nonzero edge weights
%    We pick quantiles so that heavier edges map into warmer colours.
% -------------------------------------------------------------------------
nz = nonzeros(ThisNetwork);
% Bins: 75th, 85th, 90th, 95th percentiles of nonzero weights
Qvs = [0.5 0.6 0.75 0.9];
Q = quantile(nz, Qvs);

% -------------------------------------------------------------------------
% 2) Node layout: equally spaced points on a circle of radii r = [rx ry]
% -------------------------------------------------------------------------
NumParties = size(ThisNetwork,1);
th = linspace(0, 2*pi, NumParties+1);  % last point repeats the first

% -------------------------------------------------------------------------
% 3) Optional circular background badge (white or dark)
% -------------------------------------------------------------------------
W=1;
if W == 1
    % White badge
    pp = patch(ThisYear + 1.25.*r(1).*cos(th), ...
        YS       + 1.25.*r(2).*sin(th), 'w', ...
        'EdgeColor','none');
elseif W >= 2
    % Dark blue badge
    pp = patch(ThisYear + 1.25.*r(1).*cos(th), ...
        YS       + 1.25.*r(2).*sin(th), [32 40 90]./255, ...
        'EdgeColor','none');
end

% -------------------------------------------------------------------------
% 4) Light edge pass (behind): draw thin grey connectors for small weights
%    Only if W ~= 3 (i.e., not in “nodes-only” mode).
% -------------------------------------------------------------------------
if W ~= 3
    for i = 1:NumParties
        for j = i+1:NumParties
            wij = ThisNetwork(i,j);
            if wij > 0
                % Bin index (0 means below 75th percentile -> draw as grey)
                ThisCL = sum(Q < wij);
                if ThisCL == 0
                    plot( ThisYear + r(1).*cos(th([i j])), ...
                        YS       + r(2).*sin(th([i j])), ...
                        '-', 'Color', 0.85.*ones(1,3), 'LineWidth', LW(1) );
                end
            end
        end
    end
end

% -------------------------------------------------------------------------
% 6) Heavy edge pass (front): draw coloured edges above quantile thresholds
% -------------------------------------------------------------------------
if W ~= 3
    for i = 1:NumParties
        for j = i+1:NumParties
            wij = ThisNetwork(i,j);
            if wij > 0
                ThisCL = sum(Q < wij);  % 0..4
                if ThisCL > 0
                    plot( ThisYear + r(1).*cos(th([i j])), ...
                        YS       + r(2).*sin(th([i j])), ...
                        '-', 'Color', CL(ThisCL,:), 'LineWidth', LW(2) );
                end
            end
        end
    end
end

% -------------------------------------------------------------------------
% 5) Nodes: draw party markers
% -------------------------------------------------------------------------
for i = 1:NumParties
    xi = ThisYear + r(1).*cos(th(i));
    yi = YS       + r(2).*sin(th(i));
    plot(xi, yi, '.', 'color', MKT, 'MarkerSize', MKS);
end

% -------------------------------------------------------------------------
% 7) Optional label above the glyph (with background matching the badge)
% -------------------------------------------------------------------------
if ~isempty(Title)
    tt = text(ThisYear, YS + TShift, Title, ...
        'HorizontalAlignment','center', 'FontSize', FS, 'color', MKT);
    % if W == 1
    %     set(tt, 'BackgroundColor','w');
    % elseif W >= 2
    %     set(tt, 'BackgroundColor',[32 40 90]./255, 'Color','w');
    % end
end
end