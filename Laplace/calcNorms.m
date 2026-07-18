function [L2_norm,H1_norm] = calcNorms(uSolution_cvec,connectivityMatrix_mat,connectivityMatrixGeo_mat,totNumElements,xCoordGeo_vec,yCoordGeo_vec,n_ElementType,n_GeometryElementType,order,order_G,numGeometryGaussPoints,numGaussPoints,xCoord_vec,yCoord_vec)

    L2_integral = 0;
    H1_integral = 0;

    for eleNum = 1:totNumElements

        elementNodes_vec = connectivityMatrix_mat(eleNum,:);
        xNodesVals_vec   = xCoord_vec(elementNodes_vec);
        yNodesVals_vec   = yCoord_vec(elementNodes_vec);

        elementNodesGeo_vec = connectivityMatrixGeo_mat(eleNum,:);
        xNodesValsGeo_vec   = xCoordGeo_vec(elementNodesGeo_vec);
        yNodesValsGeo_vec   = yCoordGeo_vec(elementNodesGeo_vec);

        [weights_pages,N_row_points_pages,~,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,~] = getGaussRelated(xNodesValsGeo_vec,yNodesValsGeo_vec,n_ElementType,numGaussPoints,n_GeometryElementType,numGeometryGaussPoints,order,order_G);
        
        dNdX_vec_pages = N_diff_PhysCoords_points_rows_pages(1, :, :);
        dNdY_vec_pages = N_diff_PhysCoords_points_rows_pages(2, :, :);

        % Physical coordinates at Gauss points: sum over geometry nodes
        xNodes_row = xNodesVals_vec';
        yNodes_row = yNodesVals_vec';
        
        % Compute x and y at each Gauss point
        x_gp_pages = sum(xNodes_row .* N_row_points_pages, 2);
        y_gp_pages = sum(yNodes_row .* N_row_points_pages, 2);
        
        % Evaluate exact solution at Gauss points directly
        uExact_gp  = (cosh(pi*y_gp_pages) - coth(pi)*sinh(pi*y_gp_pages)) .* sin(pi*x_gp_pages);
        uxExact_gp = pi * (cosh(pi*y_gp_pages) - coth(pi)*sinh(pi*y_gp_pages)) .* cos(pi*x_gp_pages);
        uyExact_gp = pi * (sinh(pi*y_gp_pages) - coth(pi)*cosh(pi*y_gp_pages)) .* sin(pi*x_gp_pages);

        % Evaluate FEM solution at Gauss points
        uSol_nodes = uSolution_cvec(elementNodes_vec)';

        uSol_gp  = sum(uSol_nodes .* N_row_points_pages, 2);
        uxSol_gp = sum(uSol_nodes .* dNdX_vec_pages, 2);
        uySol_gp = sum(uSol_nodes .* dNdY_vec_pages, 2);

        lError_gp = (uExact_gp - uSol_gp).^2;
        hError_gp = lError_gp + (uxExact_gp - uxSol_gp).^2 + (uyExact_gp - uySol_gp).^2;

        L2_integral = L2_integral + sum(lError_gp .* weights_pages .* J_det_points_elev, 3);
        H1_integral = H1_integral + sum(hError_gp .* weights_pages .* J_det_points_elev, 3);

    end

    L2_norm = sqrt(sum(L2_integral(:)));
    H1_norm = sqrt(sum(H1_integral(:)));
end
