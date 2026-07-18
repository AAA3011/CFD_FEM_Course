% makeVTKdir - Ensure VTK output directory exists.
%
% Description: Create output folder if it does not already exist.
% Inputs:
%   none
% Outputs:
%   none
function []=makeVTKdir(outputFolder)

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

end
