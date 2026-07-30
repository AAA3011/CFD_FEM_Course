%% ================================================================
%%  plotImplicitResults.m
%%  Professional post-processing for the implicit Newton-Raphson
%%  lid-driven cavity solver — three Reynolds numbers side by side.
%%
%%  Requires:
%%    cavity_implicit_Re10.mat
%%    cavity_implicit_Re100.mat
%%    cavity_implicit_Re1000.mat
%%
%%  Produces (saved as SVG in ./figures_implicit/):
%%    fig1_streamlines.svg          — streamlines  (1x3 panel)
%%    fig2_speed.svg                — velocity magnitude (1x3 panel)
%%    fig3_pressure.svg             — pressure (1x3 panel)
%%    fig4_u_contour.svg            — u-velocity (1x3 panel)
%%    fig5_v_contour.svg            — v-velocity (1x3 panel)
%%    fig6_u_profile.svg            — centreline u validation (1x3 panel)
%%    fig7_newton_convergence.svg   — Newton convergence per time step
%%                                    (implicit-solver specific plot)
%% ================================================================

clc; clearvars; close all;

%% ── Output folder ────────────────────────────────────────────────
outDir = 'figures_implicit';
if ~exist(outDir, 'dir'); mkdir(outDir); end

%% ── Load the three runs ──────────────────────────────────────────
cases = {10, 100, 1000};
nCases = numel(cases);
data   = cell(nCases, 1);
for k = 1:nCases
    fname = sprintf('cavity_implicit_Re%d.mat', cases{k});
    data{k} = load(fname);
    fprintf('Loaded %s\n', fname);
end

%% ── Reference data (Pakdel et al.) ──────────────────────────────
%  Each CSV is expected to have two columns: [u, y]
refFiles = { ...
    'x_velocityVariationPakdel10Re.csv', ...
    'x_velocityVariationPakdel100Re.csv', ...
    'x_velocityVariationPakdel1000Re.csv'};
refData = cell(nCases, 1);
for k = 1:nCases
    if exist(refFiles{k}, 'file')
        m = readmatrix(refFiles{k}, 'FileType','text');
        refData{k}.u = m(:,1);
        refData{k}.y = m(:,2);
    else
        refData{k} = [];
        fprintf('Reference file not found: %s — skipping markers.\n', refFiles{k});
    end
end

%% ── Shared interpolation grid ────────────────────────────────────
%% ── Shared interpolation grid ────────────────────────────────────
nPlot  = 600;   % Increased from 200 for a finer underlying grid
xPlot1 = linspace(0, 1, nPlot);
yPlot1 = linspace(0, 1, nPlot);
[Xg, Yg] = meshgrid(xPlot1, yPlot1);
%% ── Typography & colours (matching Ch-3 style) ───────────────────
set(groot, 'defaultTextInterpreter',        'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter',      'latex');
set(groot, 'defaultAxesFontSize',           11);
set(groot, 'defaultLineLineWidth',          1.4);

% FIX: all entries must be numeric RGB triplets (1x3 double). Mixing
% char color specs ('b','r') with numeric triplets caused
% [panelColors{k}; panelColors{k}] to become an invalid char array when
% k==1 or k==2, which colormap() tried (and failed) to parse as a named
% colormap string, throwing "First argument must be text." in strfind.
panelColors = {[0.00 0.00 1.00], ...   % blue
               [1.00 0.00 0.00], ...   % red
               [0.10 0.55 0.10]};      % green
reLabels    = {'Re = 10', 'Re = 100', 'Re = 1000'};

%% ── Helper: tight panel margins ─────────────────────────────────
panelW = 0.265;   panelH = 0.74;
panelX = [0.055  0.375  0.695];
panelY = 0.13;

%% ================================================================
%%  Fig 1 — Streamlines (1 × 3)
%% ================================================================
fig = figure('Units','centimeters','Position',[2 2 36 12]);
sgtitle('\textbf{Streamlines} --- Lid-Driven Cavity Flow (Implicit Solver)', ...
        'Interpreter','latex','FontSize',13);

for k = 1:nCases
    d  = data{k};
    Fu = scatteredInterpolant(d.xCoord_vec(:), d.yCoord_vec(:), d.un_col(:), ...
                              'natural','nearest');
    Fv = scatteredInterpolant(d.xCoord_vec(:), d.yCoord_vec(:), d.vn_col(:), ...
                              'natural','nearest');
    U  = Fu(Xg, Yg);
    V  = Fv(Xg, Yg);
    speed = sqrt(U.^2 + V.^2);

    ax = axes('Position', [panelX(k) panelY panelW panelH]);
    streamslice(Xg, Yg, U, V, 2.5);
    colormap(ax, [panelColors{k}; panelColors{k}]);   % monochrome lines
    set(ax.Children, 'Color', panelColors{k}, 'LineWidth', 0.8);
    axis equal tight;
    box on;
    xlabel('$x/L$','Interpreter','latex');
    if k == 1; ylabel('$y/L$','Interpreter','latex'); end
    title(reLabels{k}, 'Interpreter','latex');
    set(ax,'XLim',[0 1],'YLim',[0 1],'TickDir','in');
end
set(fig, 'Color','none', 'InvertHardcopy','off');
print(fig, fullfile(outDir,'fig1_streamlines.svg'), '-dsvg', '-painters');
fprintf('Saved fig1_streamlines.svg\n');

%% ================================================================
%%  Fig 2 — Velocity Magnitude (1 × 3)
%% ================================================================
fig = figure('Units','centimeters','Position',[2 2 36 12]);
sgtitle('\textbf{Velocity Magnitude} --- Lid-Driven Cavity Flow (Implicit Solver)', ...
        'Interpreter','latex','FontSize',13);

for k = 1:nCases
    d  = data{k};
    Fu = scatteredInterpolant(d.xCoord_vec(:), d.yCoord_vec(:), d.un_col(:), ...
                              'natural','nearest');
    Fv = scatteredInterpolant(d.xCoord_vec(:), d.yCoord_vec(:), d.vn_col(:), ...
                              'natural','nearest');
    speed = sqrt(Fu(Xg,Yg).^2 + Fv(Xg,Yg).^2);

    ax = axes('Position', [panelX(k) panelY panelW panelH]);
    contourf(Xg, Yg, speed, 80, 'LineColor','flat');
    cb = colorbar; cb.TickLabelInterpreter = 'latex';
    cb.Label.String = '$|\mathbf{u}|$'; cb.Label.Interpreter = 'latex';
    colormap(ax, parula);
    axis equal tight; box on;
    xlabel('$x/L$','Interpreter','latex');
    if k == 1; ylabel('$y/L$','Interpreter','latex'); end
    title(reLabels{k}, 'Interpreter','latex');
    set(ax,'XLim',[0 1],'YLim',[0 1],'TickDir','in');
end
set(fig, 'Color','none', 'InvertHardcopy','off');
print(fig, fullfile(outDir,'fig2_speed.svg'), '-dsvg', '-painters');
fprintf('Saved fig2_speed.svg\n');

%% ================================================================
%%  Fig 3 — Pressure (1 × 3)
%% ================================================================
fig = figure('Units','centimeters','Position',[2 2 36 12]);
sgtitle('\textbf{Pressure Distribution} --- Lid-Driven Cavity Flow (Implicit Solver)', ...
        'Interpreter','latex','FontSize',13);

for k = 1:nCases
    d  = data{k};
    Fp = scatteredInterpolant(d.xCoord_vec(:), d.yCoord_vec(:), d.pn_col(:), ...
                              'natural','nearest');
    P  = Fp(Xg, Yg);

    ax = axes('Position', [panelX(k) panelY panelW panelH]);
    contourf(Xg, Yg, P, 80, 'LineColor','flat');
    cb = colorbar; cb.TickLabelInterpreter = 'latex';
    cb.Label.String = '$p$'; cb.Label.Interpreter = 'latex';
    colormap(ax, jet);
    axis equal tight; box on;
    xlabel('$x/L$','Interpreter','latex');
    if k == 1; ylabel('$y/L$','Interpreter','latex'); end
    title(reLabels{k}, 'Interpreter','latex');
    set(ax,'XLim',[0 1],'YLim',[0 1],'TickDir','in');
end
set(fig, 'Color','none', 'InvertHardcopy','off');
print(fig, fullfile(outDir,'fig3_pressure.svg'), '-dsvg', '-painters');
fprintf('Saved fig3_pressure.svg\n');

%% ================================================================
%%  Fig 4 — u-velocity contour (1 × 3)
%% ================================================================
fig = figure('Units','centimeters','Position',[2 2 36 12]);
sgtitle('$u$\textbf{-Velocity Distribution} --- Lid-Driven Cavity Flow (Implicit Solver)', ...
        'Interpreter','latex','FontSize',13);

for k = 1:nCases
    d  = data{k};
    Fu = scatteredInterpolant(d.xCoord_vec(:), d.yCoord_vec(:), d.un_col(:), ...
                              'natural','nearest');
    U  = Fu(Xg, Yg);

    ax = axes('Position', [panelX(k) panelY panelW panelH]);
    contourf(Xg, Yg, U, 80, 'LineColor','flat');
    cb = colorbar; cb.TickLabelInterpreter = 'latex';
    cb.Label.String = '$u$'; cb.Label.Interpreter = 'latex';
    colormap(ax, coolwarm_map());
    axis equal tight; box on;
    xlabel('$x/L$','Interpreter','latex');
    if k == 1; ylabel('$y/L$','Interpreter','latex'); end
    title(reLabels{k}, 'Interpreter','latex');
    set(ax,'XLim',[0 1],'YLim',[0 1],'TickDir','in');
end
set(fig, 'Color','none', 'InvertHardcopy','off');
print(fig, fullfile(outDir,'fig4_u_contour.svg'), '-dsvg', '-painters');
fprintf('Saved fig4_u_contour.svg\n');

%% ================================================================
%%  Fig 5 — v-velocity contour (1 × 3)
%% ================================================================
fig = figure('Units','centimeters','Position',[2 2 36 12]);
sgtitle('$v$\textbf{-Velocity Distribution} --- Lid-Driven Cavity Flow (Implicit Solver)', ...
        'Interpreter','latex','FontSize',13);

for k = 1:nCases
    d  = data{k};
    Fv = scatteredInterpolant(d.xCoord_vec(:), d.yCoord_vec(:), d.vn_col(:), ...
                              'natural','nearest');
    V  = Fv(Xg, Yg);

    ax = axes('Position', [panelX(k) panelY panelW panelH]);
    contourf(Xg, Yg, V, 80, 'LineColor','flat');
    cb = colorbar; cb.TickLabelInterpreter = 'latex';
    cb.Label.String = '$v$'; cb.Label.Interpreter = 'latex';
    colormap(ax, coolwarm_map());
    axis equal tight; box on;
    xlabel('$x/L$','Interpreter','latex');
    if k == 1; ylabel('$y/L$','Interpreter','latex'); end
    title(reLabels{k}, 'Interpreter','latex');
    set(ax,'XLim',[0 1],'YLim',[0 1],'TickDir','in');
end
set(fig, 'Color','none', 'InvertHardcopy','off');
print(fig, fullfile(outDir,'fig5_v_contour.svg'), '-dsvg', '-painters');
fprintf('Saved fig5_v_contour.svg\n');

%% ================================================================
%%  Fig 6 — Centreline u-velocity validation (1 × 3)
%% ================================================================
fig = figure('Units','centimeters','Position',[2 2 36 12]);
sgtitle(['$u$\textbf{-Velocity Profile at } $x{=}0.5$ --- ' ...
         'Implicit Solver vs.\ Reference Data'], ...
        'Interpreter','latex','FontSize',13);

for k = 1:nCases
    d = data{k};

    % extract nodes on centreline x = 0.5
    tolCL = 1e-10;
    ind   = find(abs(d.xCoord_vec - 0.5) < tolCL);
    yCL   = d.yCoord_vec(ind);
    uCL   = d.un_col(ind);
    [yCL, ord] = sort(yCL);
    uCL = uCL(ord);

    ax = axes('Position', [panelX(k) panelY panelW panelH]);
    hold on; box on; grid on;
    plot(uCL, yCL, '-', 'Color', panelColors{k}, 'LineWidth', 1.8, ...
         'DisplayName','Implicit Solver');
    if ~isempty(refData{k})
        plot(refData{k}.u, refData{k}.y, 'o', ...
             'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8 0.8 0.8], ...
             'MarkerSize', 5, 'LineWidth', 0.8, ...
             'DisplayName','Pakdel et al.');
    end
    xlabel('$u$',  'Interpreter','latex');
    if k == 1; ylabel('$y/L$','Interpreter','latex'); end
    title(reLabels{k}, 'Interpreter','latex');
    legend('Location','best','FontSize',9);
    set(ax,'YLim',[0 1],'TickDir','in');
    hold off;
end
set(fig, 'Color','none', 'InvertHardcopy','off');
print(fig, fullfile(outDir,'fig6_u_profile.svg'), '-dsvg', '-painters');
fprintf('Saved fig6_u_profile.svg\n');

%% ================================================================
%%  Fig 7 — Newton Convergence History  [IMPLICIT-SOLVER SPECIFIC]
%%
%%  To use this plot you need to save the convergence history inside
%%  computeDeltaVals.m.  Add the following lines at the end of the
%%  Newton while-loop (after the norm calculation):
%%
%%      convHistory(iter) = nor;
%%
%%  and before the function ends:
%%
%%      unplus_col = uKplus_col;
%%      vnplus_col = vKplus_col;
%%      pnplus_col = pKplus_col;
%%      % already there
%%
%%  Then in main.m collect it per time step:
%%
%%      newtonHist{itr} = convHistory;   % cell array declared before loop
%%
%%  and add  'newtonHist'  to the save() call.
%%
%%  If newtonHist is not present in the .mat files the plot is skipped.
%% ================================================================
fig = figure('Units','centimeters','Position',[2 2 36 12]);
sgtitle(['\textbf{Newton--Raphson Convergence History}' ...
         ' --- Representative Time Steps'], ...
        'Interpreter','latex','FontSize',13);

anyPlotted = false;
for k = 1:nCases
    d = data{k};
    ax = axes('Position', [panelX(k) panelY panelW panelH]);
    hold on; box on; grid on;

    if isfield(d, 'newtonHist')
        nh     = d.newtonHist;
        nSteps = numel(nh);
        % pick 5 representative time steps spread across the run
        stepIdx = unique(round(linspace(1, nSteps, 5)));
        cmap    = lines(numel(stepIdx));
        for s = 1:numel(stepIdx)
            iStep = stepIdx(s);
            hist_s = nh{iStep}(:);
            semilogy(1:numel(hist_s), hist_s, '-o', ...
                'Color', cmap(s,:), 'MarkerSize', 4, 'LineWidth', 1.2, ...
                'DisplayName', sprintf('Step %d', iStep));
            anyPlotted = true;
        end
        yline(1e-12, '--k', '$\varepsilon_{\mathrm{tol}}$', ...
              'Interpreter','latex','LabelHorizontalAlignment','right', ...
              'LineWidth',1.0);
        xlabel('Newton iteration', 'Interpreter','latex');
        if k == 1; ylabel('$\|\Delta\mathbf{x}\|_2$','Interpreter','latex'); end
        title(reLabels{k}, 'Interpreter','latex');
        legend('Location','northeast','FontSize',8);
        set(ax,'TickDir','in');
    else
        text(0.5, 0.5, {'{\tt newtonHist} not found', ...
                         'See header comments', 'to enable this plot.'}, ...
             'Units','normalized','HorizontalAlignment','center', ...
             'Interpreter','latex','FontSize',10);
        axis off;
    end
    hold off;
end

set(fig, 'Color','none', 'InvertHardcopy','off');
print(fig, fullfile(outDir,'fig7_newton_convergence.svg'), '-dsvg', '-painters');
fprintf('Saved fig7_newton_convergence.svg\n');

fprintf('\nAll figures saved to ./%s/\n', outDir);

%% ================================================================
%%  LOCAL HELPER — cool-warm diverging colormap (blue→white→red)
%%  (built-in only from R2023b; this works on any release)
%% ================================================================
function cmap = coolwarm_map(n)
    if nargin < 1; n = 256; end
    top    = [0.230 0.299 0.754];   % cool blue
    mid    = [0.865 0.865 0.865];   % neutral grey-white
    bot    = [0.706 0.016 0.150];   % warm red
    half   = ceil(n/2);
    c1 = interp1([0 1], [top; mid], linspace(0,1,half));
    c2 = interp1([0 1], [mid; bot], linspace(0,1,n-half+1));
    cmap   = [c1; c2(2:end,:)];
end
