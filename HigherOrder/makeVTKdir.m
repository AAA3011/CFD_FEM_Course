function []=makeVTKdir(outputFolder)

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

end