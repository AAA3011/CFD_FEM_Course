% makeVTKdir - Ensure VTK output directory exists.
%
% FILE: makeVTKdir.m
% DESCRIPTION:
% Ensure the VTK output directory exists; create it if missing.
%
% Inputs:
%   none
% Outputs:
%   none
function []=makeVTKdir(outputFolder)

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

end
