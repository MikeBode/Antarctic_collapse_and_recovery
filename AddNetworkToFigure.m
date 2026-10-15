function AddNetworkToFigure(Network, X, r, Y, Title, TShift, FS, LW, NodeCol, MKS, CL)
% ADDNETWORKTOFIGURE  Draws one small circular network at (X, Y).
%
% Nodes sit evenly around a circle on a white disc. Edges below the median
% nonzero weight are drawn in light grey. Heavier edges are drawn on top,
% coloured by which weight quantile (50th, 60th, 75th or 90th) they exceed.
%
% Inputs
%   Network  N x N symmetric weight matrix
%   X, Y     centre of the glyph
%   r        [rx ry] radius of the node circle
%   Title    label drawn at the glyph ('' for none)
%   TShift   vertical offset of the label
%   FS       label font size
%   LW       [light heavy] edge line widths
%   NodeCol  RGB colour of the nodes and label
%   MKS      node marker size
%   CL       4 x 3 colours for the four heavy-edge bins

% Weight thresholds for the edge colour bins
Q = quantile(nonzeros(Network), [0.5 0.6 0.75 0.9]);

% Node angles around the circle
NumParties = size(Network, 1);
th = linspace(0, 2*pi, NumParties+1);

% White background disc
patch(X + 1.25.*r(1).*cos(th), Y + 1.25.*r(2).*sin(th), 'w', 'EdgeColor', 'none');

% Light edges first, so they sit behind
for i = 1:NumParties
    for j = i+1:NumParties
        wij = Network(i,j);
        if wij > 0 && sum(Q < wij) == 0
            plot(X + r(1).*cos(th([i j])), Y + r(2).*sin(th([i j])), ...
                '-', 'Color', 0.85.*ones(1,3), 'LineWidth', LW(1));
        end
    end
end

% Heavier edges, coloured by bin
for i = 1:NumParties
    for j = i+1:NumParties
        wij = Network(i,j);
        bin = sum(Q < wij);
        if wij > 0 && bin > 0
            plot(X + r(1).*cos(th([i j])), Y + r(2).*sin(th([i j])), ...
                '-', 'Color', CL(bin,:), 'LineWidth', LW(2));
        end
    end
end

% Nodes
for i = 1:NumParties
    plot(X + r(1).*cos(th(i)), Y + r(2).*sin(th(i)), '.', 'color', NodeCol, 'MarkerSize', MKS);
end

% Label
if ~isempty(Title)
    text(X, Y + TShift, Title, 'HorizontalAlignment', 'center', 'FontSize', FS, 'color', NodeCol);
end
end
