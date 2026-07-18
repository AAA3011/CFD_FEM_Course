% plot - For Preprocessing.
%
% FILE: Plot.m
% DESCRIPTION:
% Post-processing script for Re=100 cylinder flow. Loads results.mat,
% interpolates fields on grid, computes derived quantities (vorticity,
% pressure coefficient), and generates 12 SVG visualization plots.
%
% Inputs:
%   fig : MATLAB figure handle to save (e.g., gcf or a handle returned by figure()).
%   folder : Output folder path (string) where SVG files will be written.
%   name : File name or base name (string) for the saved SVG (no path required).
% Outputs:
%   none
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
