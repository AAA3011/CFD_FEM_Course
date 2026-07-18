% writeFinalVTK - Write VTK output files.
%
% FILE: writeFinalVTK.m
% DESCRIPTION:
% Convenience wrapper to write the final solution to a VTK file in the
% specified output folder using `newVTK`.
%
% Inputs:
%   none
% Outputs:
%   none
function []=writeFinalVTK(xCoord_vec,yCoord_vec,un_col,vn_col,connectivityMatrix_mat,pressureSolution_col,time,outputFolder)

    filename = sprintf('%s/flow_solution_final.vtk', outputFolder);
    newVTK(xCoord_vec, yCoord_vec, un_col, vn_col, connectivityMatrix_mat, pressureSolution_col, filename, time);
    fprintf('Final VTK file written: %s\n', filename);

end
