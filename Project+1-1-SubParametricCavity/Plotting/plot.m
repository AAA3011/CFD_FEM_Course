%   fig1_streamlines.svg   – streamlines only,        1x3 panel
%   fig2_speed.svg         – |u| contours,            1x3 panel
%   fig3_pressure.svg      – pressure,                1x3 panel
%   fig4_validation.svg     – u-profile,               1x3 panel (wide)
%   fig5_meshComparison.svg – solution mesh (row 1) vs geometry mesh
%                             (row 2), one column per Re case, 2x3 panel
%
%  Requires .mat files saved with: uSolution_col, vSolution_col,
%  pressureSolution_col, xCoord_vec, yCoord_vec, connectivityMatrix_mat,
%  xCoordGeo_vec, yCoordGeo_vec, connectivityMatrixGeo_mat
%  (geometry row shows a placeholder message if those fields aren't found)

clc; clearvars; close all;

%% ── Typography ────────────────────────────────────────────────────────────
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',13);

%% ── Settings ──────────────────────────────────────────────────────────────
ReList   = [10, 100, 1000];
matFiles = {'cavity_resultsRe10.mat', ...
            'cavity_resultsRe100.mat', ...
            'cavity_resultsRe1000.mat'};
csvFiles = {'x_velocityVariationPakdel10Re.csv', ...
            'x_velocityVariationPakdel100Re.csv', ...
            'x_velocityVariationPakdel1000Re.csv'};
Ngrid = 300;
LC    = [0.18 0.38 0.75;   % blue  Re=10
         0.80 0.15 0.15;   % red   Re=100
         0.10 0.62 0.25];  % green Re=1000

%% ── Load ──────────────────────────────────────────────────────────────────
D = struct();
for k = 1:3
    S  = load(matFiles{k});
    fn = fieldnames(S);

    xv = pf(S,fn,{'xCoord_vec','xCoord_cvec','xCoord','x_vec','x'});
    yv = pf(S,fn,{'yCoord_vec','yCoord_cvec','yCoord','y_vec','y'});
    u  = pf(S,fn,{'uSolution_col','uSolution_cvec','uSolution','u_col','u_cvec','u'});
    v  = pf(S,fn,{'vSolution_col','vSolution_cvec','vSolution','v_col','v_cvec','v'});
    p  = pf(S,fn,{'pressureSolution_col','pressureSolution_cvec','pressureSolution','p_col','p_cvec','p'});
    xv=xv(:); yv=yv(:); u=u(:); v=v(:); p=p(:);

    % Mesh connectivity (solution mesh = possibly higher-order, geometry mesh = Q4)
    % NOTE: use pfmat (no flattening) for connectivity matrices; pfopt for
    % optional fields that may be missing in older saved files.
    connSol  = pfmat(S,fn,{'connectivityMatrix_mat','connectivityMatrix','conn_mat'});
    connGeo  = pfopt(S,fn,{'connectivityMatrixGeo_mat','connectivityMatrixGeo','connGeo_mat'},true);
    xGeo     = pfopt(S,fn,{'xCoordGeo_vec','xCoordGeo'},false);
    yGeo     = pfopt(S,fn,{'yCoordGeo_vec','yCoordGeo'},false);

    D(k).xSol = xv; D(k).ySol = yv; D(k).connSol = connSol;
    D(k).xGeo = xGeo; D(k).yGeo = yGeo; D(k).connGeo = connGeo;

    [xg,yg] = meshgrid(linspace(0,1,Ngrid),linspace(0,1,Ngrid));
    ug = scatteredInterpolant(xv,yv,u,'linear','none'); ug=ug(xg,yg);
    vg = scatteredInterpolant(xv,yv,v,'linear','none'); vg=vg(xg,yg);
    pg = scatteredInterpolant(xv,yv,p,'linear','none'); pg=pg(xg,yg);

    D(k).xg=xg; D(k).yg=yg; D(k).ug=ug; D(k).vg=vg; D(k).pg=pg;
    D(k).spd = sqrt(ug.^2 + vg.^2);

    % u-profile at x=0.5
    [~,ci] = min(abs(xv-0.5));
    mask   = abs(xv-xv(ci)) < 1e-9;
    [ys,so]= sort(yv(mask)); utmp=u(mask);
    D(k).y_prof=ys; D(k).u_prof=utmp(so);

    % CSV reference
    if exist(csvFiles{k},'file')
        M=readmatrix(csvFiles{k},'FileType','text');
        D(k).ref_u=M(:,1); D(k).ref_y=M(:,2); D(k).hasRef=true;
        fprintf('Re=%-4d | CSV: %d pts\n',ReList(k),size(M,1));
    else
        D(k).hasRef=false;
        fprintf('Re=%-4d | CSV not found\n',ReList(k));
    end
end

%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 1 – Streamlines only  (1x3 panel, no |u| coloring)
%% ════════════════════════════════════════════════════════════════════════
f1  = figure('Units','centimeters','Position',[1 2 44 14]);
tl1 = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for k = 1:3
    nexttile;
    hs = streamslice(D(k).xg, D(k).yg, D(k).ug, D(k).vg, 2.5);
    set(hs,'Color',LC(k,:),'LineWidth',1.0);
    axis equal tight; xlim([0 1]); ylim([0 1]);
    xlabel('x/L','FontSize',13); ylabel('y/L','FontSize',13);
    title(sprintf('Re = %d',ReList(k)),'FontSize',14,'FontWeight','bold');
    set(gca,'Box','on'); grid on;
end
title(tl1,'Streamlines  -  Lid-Driven Cavity Flow', ...
      'FontSize',15,'FontWeight','bold','FontName','Times New Roman');
saveSVG(f1,'fig1_streamlines');

%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 2 – Speed magnitude |u|  (1x3 panel)
%% ════════════════════════════════════════════════════════════════════════
f2  = figure('Units','centimeters','Position',[1 2 44 14]);
tl2 = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for k = 1:3
    nexttile;
    contourf(D(k).xg, D(k).yg, D(k).spd, 80, 'LineColor','none');
    colormap(gca, turbo(256));
    cb = colorbar('southoutside');
    cb.Label.String='Speed  |u|'; cb.Label.FontSize=12;
    axis equal tight; xlim([0 1]); ylim([0 1]);
    xlabel('x/L','FontSize',13); ylabel('y/L','FontSize',13);
    title(sprintf('Re = %d',ReList(k)),'FontSize',14,'FontWeight','bold');
    set(gca,'Box','on'); grid on;
end
title(tl2,'Velocity Magnitude  -  Lid-Driven Cavity Flow', ...
      'FontSize',15,'FontWeight','bold','FontName','Times New Roman');
saveSVG(f2,'fig2_speed');

%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 3 – Pressure contours  (1x3 panel, no iso-lines)
%% ════════════════════════════════════════════════════════════════════════
f3  = figure('Units','centimeters','Position',[1 2 44 14]);
tl3 = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for k = 1:3
    nexttile;
    % Remove arbitrary constant, clip spikes
    pg_k = D(k).pg - mean(D(k).pg(:),'omitnan');
    sig  = std(pg_k(:),'omitnan');
    pg_k = max(min(pg_k, 3*sig), -3*sig);

    contourf(D(k).xg, D(k).yg, pg_k, 1000, 'LineColor','none');
    colormap(gca, redblue(256));
    plim = max(abs(pg_k(:)));
    if plim > 0, clim([-plim plim]); end

    cb = colorbar('southoutside');
    cb.Label.String='Pressure  p'; cb.Label.FontSize=12;
    axis equal tight; xlim([0 1]); ylim([0 1]);
    xlabel('x/L','FontSize',13); ylabel('y/L','FontSize',13);
    title(sprintf('Re = %d',ReList(k)),'FontSize',14,'FontWeight','bold');
    set(gca,'Box','on'); grid on;
end
title(tl3,'Pressure Contours  -  Lid-Driven Cavity Flow', ...
      'FontSize',15,'FontWeight','bold','FontName','Times New Roman');
saveSVG(f3,'fig3_pressure');

%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 4 – Validation: u-profile at x=0.5  (1x3 panel, wide)
%% ════════════════════════════════════════════════════════════════════════
f4  = figure('Units','centimeters','Position',[1 2 58 14]);
tl4 = tiledlayout(1,3,'TileSpacing','loose','Padding','compact');
for k = 1:3
    nexttile;
    ax = gca;
    hold(ax,'on');

    % FEM line
    line(ax, D(k).u_prof, D(k).y_prof, ...
         'Color',LC(k,:),'LineStyle','-','LineWidth',2.2, ...
         'DisplayName',sprintf('FEM  Re=%d',ReList(k)));

    % Reference markers
    if D(k).hasRef
        line(ax, D(k).ref_u, D(k).ref_y, ...
             'Color',[0.2 0.2 0.2],'LineStyle','none', ...
             'Marker','o','MarkerSize',7,'LineWidth',1.4, ...
             'MarkerFaceColor','none', ...
             'DisplayName','Reference (Pakdel)');
    end

    line(ax,[0 0],[0 1],'Color','k','LineStyle','--','LineWidth',0.8, ...
         'HandleVisibility','off');

    xlabel(ax,'u','FontSize',13);
    ylabel(ax,'y/L','FontSize',13);
    title(ax,sprintf('Re = %d',ReList(k)),'FontSize',14,'FontWeight','bold');
    legend(ax,'Location','northwest','FontSize',10);
    xlim(ax,[-0.7 1.1]); ylim(ax,[0 1]);
    set(ax,'Box','on'); grid(ax,'on');
end
title(tl4,'u-velocity Profile at x = 0.5  -  Validation', ...
      'FontSize',15,'FontWeight','bold','FontName','Times New Roman');
saveSVG(f4,'fig4_validation');

%% ════════════════════════════════════════════════════════════════════════
%  FIGURE 5 – Mesh comparison: Solution mesh (top row) vs Geometry mesh
%             (bottom row), each column = one Re case  (2x3 panel)
%% ════════════════════════════════════════════════════════════════════════
f5  = figure('Units','centimeters','Position',[1 2 44 26]);
tl5 = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% -- Row 1: Solution mesh (higher-order elements) --
for k = 1:3
    nexttile(k);
    plotMeshEdges(D(k).connSol, D(k).xSol, D(k).ySol, LC(k,:));
    axis equal tight; xlim([0 1]); ylim([0 1]);
    xlabel('x/L','FontSize',13); ylabel('y/L','FontSize',13);
    title(sprintf('Solution Mesh -- Re = %d',ReList(k)), ...
          'FontSize',13,'FontWeight','bold');
    set(gca,'Box','on'); grid on;
end

% -- Row 2: Geometry mesh (Q4 elements) --
for k = 1:3
    nexttile(3+k);
    if isempty(D(k).connGeo)
        text(0.5,0.5,'Geometry mesh not saved','HorizontalAlignment','center');
        axis([0 1 0 1]);
    else
        plotMeshEdges(D(k).connGeo, D(k).xGeo, D(k).yGeo, LC(k,:));
        axis equal tight; xlim([0 1]); ylim([0 1]);
    end
    xlabel('x/L','FontSize',13); ylabel('y/L','FontSize',13);
    title(sprintf('Geometry Mesh -- Re = %d',ReList(k)), ...
          'FontSize',13,'FontWeight','bold');
    set(gca,'Box','on'); grid on;
end

title(tl5,'Mesh Comparison  -  Solution (top) vs Geometry (bottom)', ...
      'FontSize',15,'FontWeight','bold','FontName','Times New Roman');
saveSVG(f5,'fig5_meshComparison');

disp('Done -- all figures saved as SVG.');

%% ════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
%% ════════════════════════════════════════════════════════════════════════
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
