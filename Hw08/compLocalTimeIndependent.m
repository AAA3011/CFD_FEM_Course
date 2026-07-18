function [m_lump,m_cvec,mExpantion_cvec,mShock_cvec]=compLocalTimeIndependent(xNodesVals_vec,numGaussPoints,n_ElementType)

    [gaussPoints_cvec, gaussWeights_cvec] = gaussLegendreQuad(numGaussPoints);

    coords_row_pages  = reshape(gaussPoints_cvec , 1, 1, []);
    weights_pages     = reshape(gaussWeights_cvec, 1, 1, []);

    [N_row_points_pages, N_diff_rows_points_pages] = shapeFunctions(n_ElementType, coords_row_pages);

    [J_det_points_elev, ~] = Jacobian(xNodesVals_vec, N_diff_rows_points_pages);
    N_vec_pages     = N_row_points_pages;
    N_vec_pagest    = permute(N_row_points_pages,[2 1 3]);

    m_lump  = sum( (N_vec_pagest.*N_vec_pages).* weights_pages .* J_det_points_elev,3);
    m_cvec  = sum( (N_vec_pagest).* weights_pages .* J_det_points_elev,3);
    mExpantion_cvec = sum( (N_vec_pagest).* weights_pages .* J_det_points_elev,3);
    mShock_cvec = sum( (N_vec_pagest).* weights_pages .* J_det_points_elev,3);

end