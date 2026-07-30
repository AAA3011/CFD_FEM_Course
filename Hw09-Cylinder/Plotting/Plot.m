%% PLOTCYLINDERRESULTS — Post-processing script for Re=100 cylinder flow
%
%  USAGE:  Run this script directly. It loads results.mat automatically.
%
%  REQUIRED FILE:
%    results.mat  —  must contain: Fx_col, Fy_col, un_col, vn_col,
%                    pressureSolution_col (or similar), time (scalar or vec)
%
%  Geometry files (Elements10.txt, Nodes10.txt, etc.) must be in CWD.
%  Output SVGs saved to  outputFolder  defined below.
%
%  Produces:
%    01_Mesh.svg            07_Streamlines.svg
%    02_Pressure.svg        08_Vorticity.svg
%    03_Cp.svg              09_Drag.svg
%    04_Ux.svg              10_Lift.svg
%    05_Uy.svg              11_PhasePortrait.svg
%    06_Umag.svg            12_FFT.svg
%% ── Global aesthetics ─────────────────────────────────────────────────────
% Force LaTeX interpreter globally for all text elements
close all
set(groot, 'defaultTextInterpreter',       'latex', ...
           'defaultAxesTickLabelInterpreter','latex', ...
           'defaultLegendInterpreter',     'latex', ...
           'defaultColorbarTickLabelInterpreter', 'latex', ...
           'defaultAxesFontName',          'Times New Roman', ...
           'defaultAxesFontSize',          13, ...
           'defaultLineLineWidth',         1.5, ...
           'defaultAxesLineWidth',         1.0, ...
           'defaultFigureColor',           'w');
figW = 800;  figH = 550;

%% ── Configuration ─────────────────────────────────────────────────────────
outputFolder = 'Cylinder100Re';
U_inf = 1;          % free-stream velocity
rho   = 1;          % density
D_cyl = 1;          % cylinder diameter
NGRID = 400;        % interpolation grid resolution

% Solver time-step (matches timeStep = 0.005 in main.m)
dt_sim = 0.005;

% Fraction of total run treated as initial transient (impulsive-start
% oscillations) to DISCARD from Cd/Cl statistics AND from the plotted
% window. Re=100 from rest typically needs a large chunk of the run to
% settle into clean periodic shedding — 0.20 was too early, bumped up.
transientFrac = 0.50;

%% ── Load results.mat ──────────────────────────────────────────────────────
fprintf('Loading results.mat …\n');
S = load('results.mat');

% Map variables (handles truncated / alternate pressure name)
Fx_col               = S.Fx_col;
Fy_col               = S.Fy_col;
un_col               = S.un_col;
vn_col               = S.vn_col;

% Pressure: try exact name first, then search for a field starting with 'press'
if isfield(S, 'pressureSolution_col')
    pressureSolution_col = S.pressureSolution_col;
else
    fnames = fieldnames(S);
    pmatch = fnames(startsWith(fnames, 'press'));
    if isempty(pmatch)
        error('Cannot find pressure variable in results.mat. Fields present: %s', ...
              strjoin(fnames, ', '));
    end
    pressureSolution_col = S.(pmatch{1});
    fprintf('  Using ''%s'' as pressure variable.\n', pmatch{1});
end

% Time: reconstruct if scalar or missing
if isfield(S, 'time') && numel(S.time) > 1
    time = S.time(:);
    fprintf('  Time vector loaded from file (%d steps).\n', numel(time));
else
    Nt   = numel(Fx_col);
    time = (0 : Nt-1)' * dt_sim;
    fprintf('  time was scalar/missing → reconstructed: %d steps, dt=%.4g\n', Nt, dt_sim);
end

%% ── Unpack workspace variables ────────────────────────────────────────────
Fx = Fx_col(:);
Fy = Fy_col(:);
u  = un_col(:);
v  = vn_col(:);
p  = pressureSolution_col(:);
t  = time(:);

%% ── Geometry ──────────────────────────────────────────────────────────────
fprintf('Reading geometry …\n');
[connectivityMatrix_mat, xCoord_vec, yCoord_vec, ~, ~, inlet_vec, cylinderWall_vec, ~] = readGeoAndBC();

x = xCoord_vec(:);
y = yCoord_vec(:);

%% ── Derived fields ────────────────────────────────────────────────────────
Umag = sqrt(u.^2 + v.^2);

%% ── Uniform interpolation grid ────────────────────────────────────────────
xg = linspace(min(x), max(x), NGRID);
yg = linspace(min(y), max(y), NGRID);
[XG, YG] = meshgrid(xg, yg);

Fp = scatteredInterpolant(x, y, p,    'linear', 'none');
Fu = scatteredInterpolant(x, y, u,    'linear', 'none');
Fv = scatteredInterpolant(x, y, v,    'linear', 'none');
FM = scatteredInterpolant(x, y, Umag, 'linear', 'none');

PG = Fp(XG, YG);
UG = Fu(XG, YG);
VG = Fv(XG, YG);
MG = FM(XG, YG);

% Vorticity  ωz = ∂v/∂x − ∂u/∂y
[dVdx, ~   ] = gradient(VG, xg, yg);
[~,    dUdy] = gradient(UG, xg, yg);
WZ = dVdx - dUdy;

% Reference pressure (inlet mean)
p_inf = mean(p(inlet_vec));

% Steady-state start index — skip impulsive-start transient
i0 = max(1, round(transientFrac * numel(t)));

% Force coefficients
Cd = 2 * Fx ./ (rho * U_inf^2 * D_cyl);
Cl = 2 * Fy ./ (rho * U_inf^2 * D_cyl);

% Steady-state statistics (used for printed values + mean lines)
Cd_ss   = Cd(i0:end);
Cl_ss   = Cl(i0:end);
t_ss    = t(i0:end);
Cd_mean = mean(Cd_ss);
Cl_mean = mean(Cl_ss);
Cl_amp  = (max(Cl_ss) - min(Cl_ss)) / 2;

fprintf('  Steady-state window: t = [%.2f, %.2f] (skipped first %.0f%% as transient)\n', ...
        t_ss(1), t_ss(end), transientFrac*100);
fprintf('  Mean Cd = %.4f | Mean Cl = %.4f | Cl amplitude = %.4f\n', Cd_mean, Cl_mean, Cl_amp);

% Output folder
if ~isfolder(outputFolder), mkdir(outputFolder); end

%% ══════════════════════════════════════════════════════════════════════════
%  01  MESH
%% ══════════════════════════════════════════════════════════════════════════
fprintf('[1/12] Mesh …\n');
fig = figure('Units','pixels','Position',[50 50 figW figH],'Color','w');
hold on;

Xe = xCoord_vec(connectivityMatrix_mat);
Ye = yCoord_vec(connectivityMatrix_mat);
Xn = [Xe, nan(size(Xe,1),1)]';
Yn = [Ye, nan(size(Ye,1),1)]';
patch(Xn(:), Yn(:), 'w', 'EdgeColor',[0.35 0.35 0.35], 'LineWidth',0.4);

axis equal tight;  box on;
xlabel('$x/d$','Interpreter','latex');
ylabel('$y/d$','Interpreter','latex');
title('Computational Mesh','Interpreter','latex');
saveSVG(fig, outputFolder, '01_Mesh');

%% ══════════════════════════════════════════════════════════════════════════
%  02  PRESSURE
%% ══════════════════════════════════════════════════════════════════════════
fprintf('[2/12] Pressure …\n');
fig = figure('Units','pixels','Position',[50 50 figW figH],'Color','w');
contourf(XG, YG, PG, 60, 'LineColor','none');
colormap(gca, cmap_redblue(256));
cb = colorbar;  cb.Label.String = '$p$';  cb.Label.Interpreter = 'latex';
hold on;  cylinderPatch(0.5);
axis equal tight;  box on;
xlabel('$x/d$','Interpreter','latex');
ylabel('$y/d$','Interpreter','latex');
title('Pressure Field','Interpreter','latex');
saveSVG(fig, outputFolder, '02_Pressure');

%% ══════════════════════════════════════════════════════════════════════════
%  03  PRESSURE COEFFICIENT  (with Park 1998 validation)
%% ══════════════════════════════════════════════════════════════════════════
fprintf('[3/12] Cp …\n');

% Load Park (1998) validation data
pakdel_mat  = readmatrix('PressureCoefficientPark.csv', 'FileType','text', 'Range','A1:B94');
xpark_vec   = pakdel_mat(:,1);   % angle [deg], 0-180
ypark_vec   = pakdel_mat(:,2);   % Cp

% Cylinder-wall nodes
xc = x(cylinderWall_vec);
yc = y(cylinderWall_vec);
pc = p(cylinderWall_vec);

Cp    = (pc - p_inf) ./ (0.5 * rho * U_inf^2);
theta = atan2d(yc, xc);

% Take one continuous half-arc in the mesh's native node order (this is
% what kept your original plot connected/non-jagged — no sorting/folding
% of two separate halves), then apply your sign flip + flip() to align
% with Park's 0->180 convention.
half       = floor(length(theta)/2);
theta_half = theta(half-1:end-1) * -1;
Cp_half    = flip(Cp(half-1:end-1));

fig = figure('Units','pixels','Position',[50 50 figW figH],'Color','w');
plot(theta_half, Cp_half, 'k-', 'LineWidth', 2, 'DisplayName','FEM (present, $Re=100$)');
hold on;
plot(xpark_vec, ypark_vec, 'ro', 'MarkerSize', 6, ...
    'MarkerFaceColor','r', 'DisplayName','Park (1998)');

xlabel('$\theta$ [deg]','Interpreter','latex');
ylabel('$C_p$','Interpreter','latex');
title('Pressure Coefficient Distribution --- Cylinder Surface ($Re=100$)','Interpreter','latex');
legend('Location','best','Interpreter','latex');
xlim([0 180]);  xticks(0:20:180);
grid on;  box on;

saveSVG(fig, outputFolder, '03_Cp');

%% ══════════════════════════════════════════════════════════════════════════
%  04  U-VELOCITY
%% ══════════════════════════════════════════════════════════════════════════
fprintf('[4/12] Ux …\n');
fig = figure('Units','pixels','Position',[50 50 figW figH],'Color','w');
contourf(XG, YG, UG, 60, 'LineColor','none');
colormap(gca, cmap_coolwarm(256));
cb = colorbar;  cb.Label.String = '$u/U_\infty$';  cb.Label.Interpreter='latex';
hold on;  cylinderPatch(0.5);
axis equal tight;  box on;
xlabel('$x/d$','Interpreter','latex');
ylabel('$y/d$','Interpreter','latex');
title('Streamwise Velocity $u$','Interpreter','latex');
saveSVG(fig, outputFolder, '04_Ux');

%% ══════════════════════════════════════════════════════════════════════════
%  05  V-VELOCITY
%% ══════════════════════════════════════════════════════════════════════════
fprintf('[5/12] Uy …\n');
fig = figure('Units','pixels','Position',[50 50 figW figH],'Color','w');
contourf(XG, YG, VG, 60, 'LineColor','none');
colormap(gca, cmap_coolwarm(256));
cb = colorbar;  cb.Label.String = '$v/U_\infty$';  cb.Label.Interpreter='latex';
hold on;  cylinderPatch(0.5);
axis equal tight;  box on;
xlabel('$x/d$','Interpreter','latex');
ylabel('$y/d$','Interpreter','latex');
title('Cross-Stream Velocity $v$','Interpreter','latex');
saveSVG(fig, outputFolder, '05_Uy');

%% ══════════════════════════════════════════════════════════════════════════
%  06  VELOCITY MAGNITUDE
%% ══════════════════════════════════════════════════════════════════════════
fprintf('[6/12] |U| …\n');
fig = figure('Units','pixels','Position',[50 50 figW figH],'Color','w');
contourf(XG, YG, MG, 60, 'LineColor','none');
colormap(gca, parula(256));
cb = colorbar;  cb.Label.String = '$|\mathbf{u}|/U_\infty$';  cb.Label.Interpreter='latex';
hold on;  cylinderPatch(0.5);
axis equal tight;  box on;
xlabel('$x/d$','Interpreter','latex');
ylabel('$y/d$','Interpreter','latex');
title('Velocity Magnitude','Interpreter','latex');
saveSVG(fig, outputFolder, '06_Umag');

%% ══════════════════════════════════════════════════════════════════════════
%  07  STREAMLINES (WITH ARROWS)
%% ══════════════════════════════════════════════════════════════════════════
fprintf('[7/12] Streamlines …\n');

fig = figure('Units','pixels','Position',[50 50 figW figH],'Color','w');
hold on;

% Use streamslice to generate streamlines with directional arrows.
% The density parameter controls how many lines are drawn.
hSL = streamslice(XG, YG, UG, VG, 3.5);

% Style the lines to match the blue theme from your reference image
set(hSL, 'Color', '#0072BD', 'LineWidth', 0.8);

% Overlay the cylinder
cylinderPatch(0.5);

axis equal tight;
box on;
grid on;

xlabel('$x/d$', 'Interpreter', 'latex');
ylabel('$y/d$', 'Interpreter', 'latex');
title('Streamlines', 'Interpreter', 'latex');

saveSVG(fig, outputFolder, '07_Streamlines');

%% ══════════════════════════════════════════════════════════════════════════
%  08  VORTICITY
%% ══════════════════════════════════════════════════════════════════════════
fprintf('[8/12] Vorticity …\n');
climV = max(abs(WZ(:)), [], 'omitnan') * 0.8;
fig = figure('Units','pixels','Position',[50 50 figW figH],'Color','w');
contourf(XG, YG, WZ, 80, 'LineColor','none');
clim([-climV climV]);
colormap(gca, cmap_redblue(256));
cb = colorbar;  cb.Label.String = '$\omega_z\,d/U_\infty$';  cb.Label.Interpreter='latex';
hold on;  cylinderPatch(0.5);
axis equal tight;  box on;
xlabel('$x/d$','Interpreter','latex');
ylabel('$y/d$','Interpreter','latex');
% Accented chars (á) in default LaTeX font can render as boxes/blank ---
% escaped with \'{a} so it always renders regardless of font.
title('Vorticity $\omega_z$ ','Interpreter','latex');
saveSVG(fig, outputFolder, '08_Vorticity');

%% ══════════════════════════════════════════════════════════════════════════
%  09  DRAG COEFFICIENT  (transient clipped from view)
%% ══════════════════════════════════════════════════════════════════════════


fprintf('[9/12] Cd …\n');
tSkip  = 40;   % skip first few time units -- impulsive-start spike dominates y-scale otherwise
maskCd = (t >= tSkip) & (t <= 400);
fig = figure('Units','pixels','Position',[50 50 figW figH],'Color','w');
plot(t(maskCd), Cd(maskCd), 'k-', 'LineWidth', 1.2);
hold on;
yline(Cd_mean, 'r--', ...
      sprintf('$\\bar{C}_D = %.3f$', Cd_mean), ...
      'Interpreter','latex', ...
      'LabelVerticalAlignment','top', ...
      'LabelHorizontalAlignment','left', ... % <-- Add this!
      'LineWidth',1.4);

xlabel('$t\,U_\infty/d$','Interpreter','latex');
ylabel('$C_D$','Interpreter','latex');
title('Drag Coefficient History','Interpreter','latex');
xlim([40 400]);
ylim([1.16 1.42]);
grid on;  box on;
saveSVG(fig, outputFolder, '09_Drag');

%% ══════════════════════════════════════════════════════════════════════════
%  09b  DRAG COEFFICIENT — CLOSEUP (t = 0 to 120)
%% ══════════════════════════════════════════════════════════════════════════
fprintf('[9b/12] Cd closeup …\n');
tSkip  = 40;   % skip first few time units -- impulsive-start spike dominates y-scale otherwise
maskCd = (t >= tSkip) & (t <= 200);
fig = figure('Units','pixels','Position',[50 50 figW figH],'Color','w');
plot(t(maskCd), Cd(maskCd), 'k-', 'LineWidth', 1.2);
xlabel('$t\,U_\infty/d$','Interpreter','latex');
ylabel('$C_D$','Interpreter','latex');
title('Drag Coefficient History --- Closeup','Interpreter','latex');
xlim([tSkip 200]);
ylim([1.16 1.42]);
grid on;  box on;
saveSVG(fig, outputFolder, '09b_Drag_Closeup');

%% ══════════════════════════════════════════════════════════════════════════
%  10  LIFT COEFFICIENT  (transient clipped from view)
%% ══════════════════════════════════════════════════════════════════════════
fprintf('[10/12] Cl …\n');
fig = figure('Units','pixels','Position',[50 50 figW figH],'Color','w');
plot(t, Cl, 'k-', 'LineWidth', 1.2);
xlabel('$t\,U_\infty/d$','Interpreter','latex');
ylabel('$C_L$','Interpreter','latex');
title('Lift Coefficient History','Interpreter','latex');
xlim([t(1) t(end)]);
grid on;  box on;
saveSVG(fig, outputFolder, '10_Lift');

%% ══════════════════════════════════════════════════════════════════════════
%  10b  LIFT COEFFICIENT — CLOSEUP (t = 0 to 120)
%% ══════════════════════════════════════════════════════════════════════════
fprintf('[10b/12] Cl closeup …\n');
tSkip = 50;
maskCl = (t >= tSkip) & (t <= 200);
fig = figure('Units','pixels','Position',[50 50 figW figH],'Color','w');
plot(t(maskCl), Cl(maskCl), 'k-', 'LineWidth', 1.2);
xlabel('$t\,U_\infty/d$','Interpreter','latex');
ylabel('$C_L$','Interpreter','latex');
title('Lift Coefficient History --- Closeup','Interpreter','latex');
xlim([tSkip 100]);
ylim([-0.32 0.32]);
grid on;  box on;
saveSVG(fig, outputFolder, '10b_Lift_Closeup');

%% ══════════════════════════════════════════════════════════════════════════
%  HELPER FUNCTIONS
%% ══════════════════════════════════════════════════════════════════════════

function saveSVG(fig, folder, name)
    fpath = fullfile(folder, [name '.svg']);
    print(fig, fpath, '-dsvg', '-r0');
    fprintf('  Saved %s\n', fpath);
end

function cylinderPatch(R)
    th = linspace(0, 2*pi, 360);
    fill(R*cos(th), R*sin(th), [0.75 0.75 0.75], ...
         'EdgeColor','k', 'LineWidth',1.0, 'HandleVisibility','off');
end

function C = cmap_redblue(n)
    if nargin < 1, n = 256; end
    h  = floor(n/2);
    lo = [linspace(0.0,1,h)', linspace(0.3,1,h)', ones(h,1)      ];
    hi = [ones(h,1),           linspace(1,0.2,h)', linspace(1,0,h)'];
    C  = [lo; hi];
end

function C = cmap_coolwarm(n)
    if nargin < 1, n = 256; end
    h    = floor(n/2);
    cool = [linspace(0.02,1,h)', linspace(0.44,1,h)', linspace(0.69,1,h)'];
    warm = [linspace(1,0.70,h)', linspace(1,0.09,h)', linspace(1,0.13,h)'];
    C    = [cool; warm];
end
