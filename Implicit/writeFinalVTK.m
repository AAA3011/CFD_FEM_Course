function []=writeFinalVTK(xCoord_vec,yCoord_vec,un_cvec,vn_cvec,connectivityMatrix_mat,pressureSolution_cvec,time,outputFolder)

    filename = sprintf('%s/flow_solution_final.vtk', outputFolder);
    newVTK(xCoord_vec, yCoord_vec, un_cvec, vn_cvec, connectivityMatrix_mat, pressureSolution_cvec, filename, time);
    fprintf('Final VTK file written: %s\n', filename);

end