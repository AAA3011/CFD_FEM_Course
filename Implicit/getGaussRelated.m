function [weights_pages,N_row_points_pages,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,ElementPhysCoords_mat] = getGaussRelated(xNodesVals_vec,yNodesVals_vec,K,numGaussPoints,K_G,numGeometryGaussPoints,order,order_G)
    %% Solution Field Part
    [gaussPoints_cvec, gaussWeights_cvec] = gaussLegendreQuad(numGaussPoints);
    
    [xi_mat, eta_mat] = meshgrid(gaussPoints_cvec, gaussPoints_cvec);    % make all different combinations of xi , eta
    coords_row_pages  = reshape([xi_mat(:), eta_mat(:)]', 1, 2, []);     % So now we have the coordinates for all compinations stored page wise
    weights_mat       = gaussWeights_cvec(:) * gaussWeights_cvec(:)';    % Make a weights tensor which means multibled weights for each coordinates set
    weights_row       = weights_mat(:)';                                 % So now we have the cooresponding weights stored in a row
    weights_pages     = reshape(weights_row,1,1,[]);

    xNodesVals_vec = xNodesVals_vec(:);
    yNodesVals_vec = yNodesVals_vec(:);
    ElementPhysCoords_mat = [xNodesVals_vec'; yNodesVals_vec'];
    
    %% Geometry Field Part
    [gaussPointsGeo_cvec, gaussWeightsGeo_cvec] = gaussLegendreQuad(numGeometryGaussPoints);
    
    [xiGeo_mat, etaGeo_mat] = meshgrid(gaussPointsGeo_cvec, gaussPointsGeo_cvec);
    coordsGeo_row_pages     = reshape([xiGeo_mat(:), etaGeo_mat(:)]', 1, 2, []);
    weightsGeo_mat          = gaussWeightsGeo_cvec(:) * gaussWeightsGeo_cvec(:)';
    weightsGeo_row          = weightsGeo_mat(:)';
    weightsGeo_pages        = reshape(weightsGeo_row,1,1,[]);
    
    %% Solution field shape function
    [ind_i, ind_j, ~] = getLagrangeIndices(order);
    % showLagrangeOrdering(ind_i, ind_j, all_nodes, order);
    
    N_row_points_pages       = N_row2D(K,ind_i,ind_j,coords_row_pages);
    % For Jacobian and physical gradients we need [dN/dxi; dN/deta].
    N_diff_rows_points_pages = N_diff_rows2D(K,ind_i,ind_j,coords_row_pages);
    
    %% Geometry shape function
    [ind_i_G, ind_j_G, ~] = getLagrangeIndices(order_G);
    % showLagrangeOrdering(ind_i_G, ind_j_G, all_nodesG, order_G);
    
    N_diff_Geo_rows_points_pages  = N_diff_rows2D(K_G,ind_i_G,ind_j_G,coordsGeo_row_pages);
    
    %% Jacobian calculation and physical derivatives
    [J_det_points_elev, N_diff_PhysCoords_points_rows_pages] = Jacobian(ElementPhysCoords_mat, N_diff_rows_points_pages,N_diff_Geo_rows_points_pages);
    % area     = sum(J_det_points_elev(:) .* weights_pages(:));
    % area_geo = sum(J_det_points_elev(:) .* weightsGeo_pages(:));

end