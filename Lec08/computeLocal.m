function [k_mat,m_mat] = computeLocal(xNodesVals_vec)
    
    numGaussPoints = 5;
    [gaussPoints_cvec, gaussWeights_cvec] = gaussLegendreQuad(numGaussPoints);

    coords_row_pages  = reshape(gaussPoints_cvec, 1, 1, []);
    weights_pages     = reshape(gaussWeights_cvec,1, 1, []);

    n_ElementType = 0;
    [N_row_points_pages, N_diff_rows_points_pages] = shapeFunctions(n_ElementType, coords_row_pages);

    [J_det_points_elev, N_diff_PhysCoords_points_rows_pages] = Jacobian(xNodesVals_vec, N_diff_rows_points_pages);

    dNdX_vec_pages = N_diff_PhysCoords_points_rows_pages;

    k_mat    = sum(permute(dNdX_vec_pages,[2 1 3])     .*  dNdX_vec_pages     .* weights_pages .* J_det_points_elev,3);
    m_mat    = sum(permute(N_row_points_pages,[2 1 3]) .*  N_row_points_pages .* weights_pages .* J_det_points_elev,3);

end
