% makeVTKdir - Ensure VTK output directory exists.
%
%
% Creates the directory if it does not exist.
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
