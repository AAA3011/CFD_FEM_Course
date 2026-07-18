% calcLocal - Calc Local.
%
% FILE: calcLocal.m
% DESCRIPTION:
% Compute the local element stiffness matrix (Laplace operator) by
% integrating derivatives of shape functions over Gauss points.
%
% Inputs:
%   elementData (variable): Struct with precomputed element integration data (weights, shape functions, derivatives, physCoords)
% Outputs:
%   none
function pk_mat = calcLocal(elementData)

    % [weights_pages,~,~,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,~] = getGaussRelated_SubParametric(xNodesValsGeo_vec,yNodesValsGeo_vec,n_ElementType,numGaussPoints,n_GeometryElementType,numGeometryGaussPoints,order,order_G);
    
    dNdX_vec_pages    = elementData.dNdX;
    dNdY_vec_pages    = elementData.dNdY;
    weights_pages     = elementData.weights;
    J_det_points_elev = elementData.J_det;
    
    % dNdX_vec_pages = N_diff_PhysCoords_points_rows_pages(1, :, :);
    % dNdY_vec_pages = N_diff_PhysCoords_points_rows_pages(2, :, :);
    
    pk_mat      = sum((pagemtimes(permute(dNdX_vec_pages,[2 1 3]),dNdX_vec_pages) + pagemtimes(permute(dNdY_vec_pages,[2 1 3]),dNdY_vec_pages))  .* weights_pages .* J_det_points_elev,3);

end
