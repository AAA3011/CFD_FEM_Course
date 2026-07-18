function [epsiK_mat,velK_mat,uRHSMult_mat,vRHSMult_mat] = computeLocalStiffMatrix(xNodesVals_vec,yNodesVals_vec,c)
    
    numGaussPoints = 4;
    [gaussPoints_cvec, gaussWeights_cvec] = gaussLegendreQuad(numGaussPoints);
    
    [xi_mat, eta_mat] = meshgrid(gaussPoints_cvec, gaussPoints_cvec);    % make all different combinations of xi , eta
    coords_row_pages  = reshape([xi_mat(:), eta_mat(:)]', 1, 2, []);     % So now we have the coordinates for all compinations stored page wise
    weights_mat       = gaussWeights_cvec(:) * gaussWeights_cvec(:)';    % Make a weights tensor which means multibled weights for each coordinates set
    weights_row       = weights_mat(:)';                                 % So no we have the cooresponding weights stored in a row
    weights_pages     = reshape(weights_row,1,1,[]);
    
    % Compute natural shape functions and their derivatives
    n_ElementType = 1;
    ElementPhysCoords_mat = [xNodesVals_vec; yNodesVals_vec]';
    [N_row_points_pages, N_diff_rows_points_pages] = shapeFunctions(n_ElementType, coords_row_pages);
    
    % Compute Jacobian and physical derivatives
    [J_det_points_elev, N_diff_PhysCoords_points_rows_pages] = Jacobian(ElementPhysCoords_mat, N_diff_rows_points_pages);
    
    
    dNdX_vec_pages = N_diff_PhysCoords_points_rows_pages(1, :, :);
    dNdY_vec_pages = N_diff_PhysCoords_points_rows_pages(2, :, :);
    
    epsiK_mat    = sum((pagemtimes(permute(N_diff_PhysCoords_points_rows_pages,[2 1 3]),N_diff_PhysCoords_points_rows_pages) + c * pagemtimes(permute(N_row_points_pages,[2 1 3]),N_row_points_pages))  .* weights_pages .* J_det_points_elev,3);
    velK_mat     = sum(pagemtimes(permute(N_row_points_pages,[2 1 3]),N_row_points_pages)  .* weights_pages .* J_det_points_elev,3);

    uRHSMult_mat    = sum(pagemtimes(permute(N_row_points_pages,[2 1 3]),dNdY_vec_pages)  .* weights_pages .* J_det_points_elev,3);
    vRHSMult_mat    = sum(pagemtimes(permute(N_row_points_pages,[2 1 3]),dNdX_vec_pages)  .* weights_pages .* J_det_points_elev,3);

end
