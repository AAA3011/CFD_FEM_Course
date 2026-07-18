% writeVTKatInterval - Write VTK output files.
%
% FILE: writeVTKatInterval.m
% DESCRIPTION:
% Write VTK output for a single time/iteration index into the specified
% output folder (used to export snapshots at intervals during simulation).
%
% Inputs:
%   none
% Outputs:
%   none
function []=writeVTKatInterval(outputFolder,itr,xCoord_vec,yCoord_vec,un_col,vn_col,connectivityMatrix_mat,pressureSolution_col,time)

    filename = sprintf('%s/flow_solution_%06d.vtk', outputFolder, itr);
    newVTK(xCoord_vec, yCoord_vec, un_col, vn_col, connectivityMatrix_mat, pressureSolution_col, filename, time);
    fprintf('Written VTK file: %s at time = %.4f, iteration = %d\n', filename, time, itr);

end
