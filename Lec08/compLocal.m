function [r_cvec,m_mat,rShock_cvec,rExpantion_cvec] = compLocal(xNodesVals_vec,un_cvec,mu,timeStep,elementLength,uShock_cvec,uExpantion_cvec)
    
    numGaussPoints = 4;
    [gaussPoints_cvec, gaussWeights_cvec] = gaussLegendreQuad(numGaussPoints);

    coords_row_pages  = reshape(gaussPoints_cvec , 1, 1, []);
    weights_pages     = reshape(gaussWeights_cvec, 1, 1, []);

    n_ElementType = 0;
    [N_row_points_pages, N_diff_rows_points_pages] = shapeFunctions(n_ElementType, coords_row_pages);

    [J_det_points_elev, N_diff_PhysCoords_points_rows_pages] = Jacobian(xNodesVals_vec, N_diff_rows_points_pages);
    dNdX_vec_pages  = N_diff_PhysCoords_points_rows_pages;
    N_vec_pages     = N_row_points_pages;

    numPages        = size(N_diff_PhysCoords_points_rows_pages,3);

    %% Initial Condtions in pages
    un_cvec         = repmat(un_cvec        ,1,1,numPages);
    uShock_cvec     = repmat(uShock_cvec    ,1,1,numPages);
    uExpantion_cvec = repmat(uExpantion_cvec,1,1,numPages);

    N_vec_pagest    = permute(N_row_points_pages,[2 1 3]);

    %% Normal Case
    u       = pagemtimes(N_vec_pages,un_cvec);
    ux      = pagemtimes(dNdX_vec_pages,un_cvec);
    tu      = 1./sqrt((2/timeStep)^2 + (2*u/elementLength).^2 + (4*mu/(elementLength^2))^2);

    W       = N_vec_pages    + tu.*u.*dNdX_vec_pages;
    wx      = dNdX_vec_pages + tu.*dNdX_vec_pages.*ux;
    
    r_cvec  = sum( (ux.*(u.*W+mu.*wx)).* weights_pages .* J_det_points_elev,3)';
    m_mat   = sum( (W.*N_vec_pagest).* weights_pages .* J_det_points_elev,3);

    %% Shock Case
    uShock      = pagemtimes(N_vec_pages,uShock_cvec);
    uxShock     = pagemtimes(dNdX_vec_pages,uShock_cvec);
    tuShock     = 1./sqrt((2/timeStep)^2 + (2*uShock/elementLength).^2 + (4*mu/(elementLength^2))^2);

    WShock      = N_vec_pages    + tuShock.*uShock.*dNdX_vec_pages;
    wxShock     = dNdX_vec_pages + tuShock.*dNdX_vec_pages.*uxShock;
    
    rShock_cvec = sum( (uxShock.*(uShock.*WShock+mu.*wxShock)).* weights_pages .* J_det_points_elev,3)';
    %% Expantion Case
    uExpantion      = pagemtimes(N_vec_pages,uExpantion_cvec);
    uxExpantion     = pagemtimes(dNdX_vec_pages,uExpantion_cvec);
    tuExpantion     = 1./sqrt((2/timeStep)^2 + (2*uExpantion/elementLength).^2 + (4*mu/(elementLength^2))^2);

    WExpantion      = N_vec_pages    + tuExpantion.*uExpantion.*dNdX_vec_pages;
    wxExpantion     = dNdX_vec_pages + tuExpantion.*dNdX_vec_pages.*uxExpantion;
    
    rExpantion_cvec = sum( (uxExpantion.*(uExpantion.*WExpantion+mu.*wxExpantion)).* weights_pages .* J_det_points_elev,3)';

end
