% plot - For Preprocessing.
%
% fig2_speed.svg         – |u| contours,            1x3 panel
% fig3_pressure.svg      – pressure,                1x3 panel
% fig4_validation.svg     – u-profile,               1x3 panel (wide)
% fig5_meshComparison.svg – solution mesh (row 1) vs geometry mesh
% (row 2), one column per Re case, 2x3 panel
%
% Requires .mat files saved with: uSolution_col, vSolution_col,
% pressureSolution_col, xCoord_vec, yCoord_vec, connectivityMatrix_mat,
% xCoordGeo_vec, yCoordGeo_vec, connectivityMatrixGeo_mat
% (geometry row shows a placeholder message if those fields aren't found)
% plot.m
% Description: Create and save validation and visualization figures from
% previously saved simulation .mat files (u, v, p, meshes).
% Inputs:
%   fig : MATLAB figure handle to save (e.g., gcf or a handle returned by figure()).
%   name : File name or base name (string) for the saved SVG (no path required).
% Outputs:
%   none
function saveSVG(fig, name)
    % Uses print for reliable SVG output
    print(fig, name, '-dsvg', '-r0');
    fprintf('Saved  %s.svg\n', name);
end

function val = pf(S, fn, candidates)
    fl = cellfun(@lower, fn, 'UniformOutput', false);
    for c = 1:numel(candidates)
        idx = find(strcmp(fl, lower(candidates{c})),1);
        if ~isempty(idx), val = S.(fn{idx})(:); return; end
    end
    error('Field not found. Tried: %s', strjoin(candidates,', '));
end

function val = pfmat(S, fn, candidates)
    % Like pf, but does NOT flatten -- use for connectivity matrices
    % (numElements x nodesPerElement), where shape must be preserved.
    fl = cellfun(@lower, fn, 'UniformOutput', false);
    for c = 1:numel(candidates)
        idx = find(strcmp(fl, lower(candidates{c})),1);
        if ~isempty(idx), val = S.(fn{idx}); return; end
    end
    error('Field not found. Tried: %s', strjoin(candidates,', '));
end

function val = pfopt(S, fn, candidates, isMatrix)
    % Like pf/pfmat, but returns [] instead of erroring if the field is
    % missing (for optional fields like geometry mesh data that may not
    % exist in older saved .mat files).
    fl = cellfun(@lower, fn, 'UniformOutput', false);
    for c = 1:numel(candidates)
        idx = find(strcmp(fl, lower(candidates{c})),1);
        if ~isempty(idx)
            if isMatrix, val = S.(fn{idx}); else, val = S.(fn{idx})(:); end
            return;
        end
    end
    val = [];
end

function plotMeshEdges(conn, x, y, color)
    % Draws element outlines for a structured quad mesh (works for both
    % Q4 geometry elements and higher-order solution elements). Node
    % ordering within an element isn't assumed -- the outline is drawn
    % via the convex hull of each element's nodes, which gives the
    % correct quad boundary for both straight-edged and higher-order
    % (curved-looking) elements on a structured mesh.
    hold on;
    numEl = size(conn,1);
    for e = 1:numEl
        nodes = conn(e,:);
        xe = x(nodes); ye = y(nodes);
        if numel(xe) >= 3
            k = convhull(xe(:), ye(:));
            patch(xe(k), ye(k), 'w', 'FaceColor','none', ...
                  'EdgeColor', color, 'LineWidth', 0.8);
        end
    end
    line(x, y, 'LineStyle','none', 'Marker','.', ...
         'Color', color*0.6, 'MarkerSize', 4);
    hold off;
end

function cmap = redblue(n)
    if nargin<1, n=256; end
    h  = floor(n/2);
    b  = [linspace(0.17,1,h)', linspace(0.27,1,h)', linspace(0.73,1,h)'];
    r  = [linspace(1,0.70,n-h)', linspace(1,0.02,n-h)', linspace(1,0.12,n-h)'];
    cmap = [b; r];
end
