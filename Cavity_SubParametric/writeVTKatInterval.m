function []=writeVTKatInterval(outputFolder,itr,xCoord_vec,yCoord_vec,un_cvec,vn_cvec,connectivityMatrix_mat,pressureSolution_cvec,time)

    filename = sprintf('%s/flow_solution_%06d.vtk', outputFolder, itr);
    newVTK(xCoord_vec, yCoord_vec, un_cvec, vn_cvec, connectivityMatrix_mat, pressureSolution_cvec, filename, time);
    fprintf('Written VTK file: %s at time = %.4f, iteration = %d\n', filename, time, itr);

end