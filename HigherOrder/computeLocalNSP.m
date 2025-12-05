function [pk_mat,umlv_cvec,vmlv_cvec] = computeLocalNSP(numGaussPoints,n_ElementType,~,~,n_GeometryElementType,numGeometryGaussPoints,xNodesValsGeo_vec,yNodesValsGeo_vec)

    [weights_pages,N_row_points_pages,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,~] = getGaussRelated(xNodesValsGeo_vec,yNodesValsGeo_vec,n_ElementType,numGaussPoints,n_GeometryElementType,numGeometryGaussPoints); 
    dNdX_vec_pages = N_diff_PhysCoords_points_rows_pages(1, :, :);
    dNdY_vec_pages = N_diff_PhysCoords_points_rows_pages(2, :, :);

    pk_mat      = sum((pagemtimes(permute(dNdX_vec_pages,[2 1 3]),dNdX_vec_pages) + pagemtimes(permute(dNdY_vec_pages,[2 1 3]),dNdY_vec_pages))  .* weights_pages .* J_det_points_elev,3);
    umlv_cvec   = sum(N_row_points_pages .* weights_pages .* J_det_points_elev,3);
    vmlv_cvec   = sum(N_row_points_pages .* weights_pages .* J_det_points_elev,3);
end