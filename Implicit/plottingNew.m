% plottingNew - Plotting New.
%
% FILE: plottingNew.m
% DESCRIPTION:
% Post-processing and plotting utilities for simulation results: loads
% saved MAT results and reference CSVs, and produces publication-ready
% figures (streamlines, contours, profiles, convergence plots).
%
% Inputs:
%   n : Number of discrete utilities or default argument placeholder (see specific functions below).
% Outputs:
%   none
function cmap = coolwarm_map(n)
% Syntax: [cmap] = coolwarm_map(n)
%
% Inputs:
%   n : Number of discrete colors to generate (integer, default 256).
%
% Outputs:
%   cmap : n-by-3 RGB colormap array.
    if nargin < 1; n = 256; end
    top    = [0.230 0.299 0.754];   % cool blue
    mid    = [0.865 0.865 0.865];   % neutral grey-white
    bot    = [0.706 0.016 0.150];   % warm red
    half   = ceil(n/2);
    c1 = interp1([0 1], [top; mid], linspace(0,1,half));
    c2 = interp1([0 1], [mid; bot], linspace(0,1,n-half+1));
    cmap   = [c1; c2(2:end,:)];
end
