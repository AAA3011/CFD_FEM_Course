% writeFinalVTK - Write VTK output files.
%
% Description: Helper to write the final simulation result VTK file.
% Inputs:
%   none
% Outputs:
%   none
function []=writeFinalVTK(xCoord_vec,yCoord_vec,un_col,vn_col,connectivityMatrix_mat,pressureSolution_col,time,outputFolder)

    filename = sprintf('%s/flow_solution_final.vtk', outputFolder);
    newVTK(xCoord_vec, yCoord_vec, un_col, vn_col, connectivityMatrix_mat, pressureSolution_col, filename, time);
    fprintf('Final VTK file written: %s\n', filename);

end
