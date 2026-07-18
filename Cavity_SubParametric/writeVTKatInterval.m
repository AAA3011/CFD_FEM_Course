% writeVTKatInterval - Write VTK output files.
%
% Description: Write a VTK file for the current iteration/time step.
% Inputs:
%   none
% Outputs:
%   none
function []=writeVTKatInterval(outputFolder,itr,xCoord_vec,yCoord_vec,un_col,vn_col,connectivityMatrix_mat,pressureSolution_col,time)

    filename = sprintf('%s/flow_solution_%06d.vtk', outputFolder, itr);
    newVTK(xCoord_vec, yCoord_vec, un_col, vn_col, connectivityMatrix_mat, pressureSolution_col, filename, time);
    fprintf('Written VTK file: %s at time = %.4f, iteration = %d\n', filename, time, itr);

end
