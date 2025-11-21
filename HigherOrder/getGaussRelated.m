function [weights_pages,N_row_points_pages,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,ElementPhysCoords_mat] = getGaussRelated(xNodesVals_vec,yNodesVals_vec,n_ElementType,numGaussPoints) 
    
    [gaussPoints_cvec, gaussWeights_cvec] = gaussLegendreQuad(numGaussPoints);
    
    [xi_mat, eta_mat] = meshgrid(gaussPoints_cvec, gaussPoints_cvec);    % make all different combinations of xi , eta
    coords_row_pages  = reshape([xi_mat(:), eta_mat(:)]', 1, 2, []);     % So now we have the coordinates for all compinations stored page wise
    weights_mat       = gaussWeights_cvec(:) * gaussWeights_cvec(:)';    % Make a weights tensor which means multibled weights for each coordinates set
    weights_row       = weights_mat(:)';                                 % So no we have the cooresponding weights stored in a row
    weights_pages     = reshape(weights_row,1,1,[]);
    
    % Compute natural shape functions and their derivatives
    ElementPhysCoords_mat = [xNodesVals_vec; yNodesVals_vec];
    ElementPhysCoords_mat = reshape(ElementPhysCoords_mat,2,length(xNodesVals_vec));

    [N_row_points_pages, N_diff_rows_points_pages] = shapeFunctions(n_ElementType, coords_row_pages);
    % N_row_points_pages = permute(N_row_points_pages,[2 1 3]);
    % Compute Jacobian and physical derivatives
    [J_det_points_elev, N_diff_PhysCoords_points_rows_pages] = Jacobian(ElementPhysCoords_mat, N_diff_rows_points_pages);
    % N_diff_PhysCoords_points_rows_pages = 
    % area = sum(J_det_points_elev(:) .* weights_pages(:))

end