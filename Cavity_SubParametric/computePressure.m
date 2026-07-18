function[pressureSolution_cvec] = computePressure(numGaussPoints,n_ElementType,connectivityMatrix_mat,xCoord_vec,yCoord_vec,pBCNodes_vec,pK_mat,totNumNodes,un_cvec,vn_cvec,epslon,n_GeometryElementType,numGeometryGaussPoints,connectivityMatrixGeo_mat)
    %% Gauss Points ans Weights Preparation
    % numGaussPoints = 4;
    % n_ElementType = 1;
    %% pRHS initialization
    pRHS_cvec = zeros(totNumNodes,1);
    for elementNumber = 1:size(connectivityMatrix_mat,1)

        elementNodesGeo_vec = connectivityMatrixGeo_mat(elementNumber,:);
        xNodesValsGeo_vec   = xCoord_vec(elementNodesGeo_vec);
        yNodesValsGeo_vec   = yCoord_vec(elementNodesGeo_vec);

        elementNodes_vec = connectivityMatrix_mat(elementNumber,:);
        % xNodesVals_vec   = xCoord_vec(elementNodes_vec);
        % yNodesVals_vec   = yCoord_vec(elementNodes_vec);
        unNodes_cvec     = un_cvec(elementNodes_vec,1);
        vnNodes_cvec     = vn_cvec(elementNodes_vec,1);
    
        [weights_pages,N_row_points_pages,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,~] = getGaussRelated(xNodesValsGeo_vec,yNodesValsGeo_vec,n_ElementType,numGaussPoints,n_GeometryElementType,numGeometryGaussPoints);
        dNdX_vec_pages = N_diff_PhysCoords_points_rows_pages(1, :, :);
        dNdY_vec_pages = N_diff_PhysCoords_points_rows_pages(2, :, :);
    
        unx_pages = pagemtimes(dNdX_vec_pages,unNodes_cvec);
        vny_pages = pagemtimes(dNdY_vec_pages,vnNodes_cvec);
    
        prhs_cvec  = sum(((-1/epslon).*(unx_pages+vny_pages).*permute(N_row_points_pages,[2 1 3]))  .* weights_pages .* J_det_points_elev,3);
        
        pRHS_cvec(elementNodes_vec) = pRHS_cvec(elementNodes_vec) + prhs_cvec;
    end
    
    %% Pressure Boundary Conditions
    pRHS_cvec(pBCNodes_vec)  = 0;
    
    %% Pressure Solution
    pressureSolution_cvec = pK_mat \ pRHS_cvec;
end