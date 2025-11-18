function [urhs_cvec,vrhs_cvec,ReLocal,cfl,xCenter,yCenter,uxrhs_cvec,uyrhs_cvec,vxrhs_cvec,vyrhs_cvec] = computeLocalNS(numGaussPoints,n_ElementType,xNodesVals_vec,yNodesVals_vec,unNodes_cvec,vnNodes_cvec,pressureSolution_cvec,Re,timeStep)
 
    % numGaussPoints = 4;
    % n_ElementType  = 1;
    [weights_pages,N_row_points_pages,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,ElementPhysCoords_mat] = getGaussRelated(xNodesVals_vec,yNodesVals_vec,n_ElementType,numGaussPoints);
    
    xphysCoord_vec = ElementPhysCoords_mat(:,1)';
    yphysCoord_vec = ElementPhysCoords_mat(:,2)';

    dNdX_vec_pages = N_diff_PhysCoords_points_rows_pages(1, :, :);
    dNdY_vec_pages = N_diff_PhysCoords_points_rows_pages(2, :, :);
 
    unx_pages  = pagemtimes(dNdX_vec_pages,unNodes_cvec);
    uny_pages  = pagemtimes(dNdY_vec_pages,unNodes_cvec);
    vny_pages  = pagemtimes(dNdY_vec_pages,vnNodes_cvec);
    vnx_pages  = pagemtimes(dNdX_vec_pages,vnNodes_cvec);
    un_pages   = pagemtimes(N_row_points_pages,unNodes_cvec);
    vn_pages   = pagemtimes(N_row_points_pages,vnNodes_cvec);
    pnx_pages  = pagemtimes(dNdX_vec_pages,pressureSolution_cvec);
    pny_pages  = pagemtimes(dNdY_vec_pages,pressureSolution_cvec);

    le1            = sqrt((yphysCoord_vec(3)-yphysCoord_vec(1))^2+(xphysCoord_vec(3)-xphysCoord_vec(1))^2);
    le2            = sqrt((yphysCoord_vec(4)-yphysCoord_vec(2))^2+(xphysCoord_vec(4)-xphysCoord_vec(2))^2);
    elementLength  = max(le1,le2);

    Ve_pages      = sqrt(un_pages.^2 + vn_pages.^2);
    thu           = 1./sqrt((2/(timeStep))^2 + ((2*Ve_pages)/elementLength).^2 + (4./((elementLength^2)*Re)^2));
    W_cvec_pages  = permute(N_row_points_pages,[2 1 3]) + thu.*(un_pages.*(permute(dNdX_vec_pages,[2 1 3])) + vn_pages.*(permute(dNdY_vec_pages,[2 1 3])));

    u             = sum(unNodes_cvec)/4;
    v             = sum(vnNodes_cvec)/4;

    uCenter      = sqrt((u)^2 + (v)^2);
    ReLocal      = uCenter * Re * elementLength;
    cfl          = uCenter * timeStep / elementLength;
    xCenter      = sum(xphysCoord_vec)/4;
    yCenter      = sum(yphysCoord_vec)/4;
    urhs_cvec    = sum(((unNodes_cvec.*permute(N_row_points_pages,[2 1 3]))- (timeStep * (un_pages.*unx_pages + vn_pages.*uny_pages) .* W_cvec_pages) - (timeStep * pnx_pages .* permute(N_row_points_pages,[2 1 3])) - ((timeStep/Re)*(unx_pages.*(permute(dNdX_vec_pages,[2 1 3])) + uny_pages.*(permute(dNdY_vec_pages,[2 1 3]))))) .* weights_pages .* J_det_points_elev,3);
    vrhs_cvec    = sum(((vnNodes_cvec.*permute(N_row_points_pages,[2 1 3]))- (timeStep * (un_pages.*vnx_pages + vn_pages.*vny_pages) .* W_cvec_pages) - (timeStep * pny_pages .* permute(N_row_points_pages,[2 1 3])) - ((timeStep/Re)*(vnx_pages.*(permute(dNdX_vec_pages,[2 1 3])) + vny_pages.*(permute(dNdY_vec_pages,[2 1 3]))))) .* weights_pages .* J_det_points_elev,3);
    uxrhs_cvec   = sum((N_row_points_pages.*unx_pages).* weights_pages .* J_det_points_elev,3);
    uyrhs_cvec   = sum((N_row_points_pages.*uny_pages).* weights_pages .* J_det_points_elev,3);
    vxrhs_cvec   = sum((N_row_points_pages.*vnx_pages).* weights_pages .* J_det_points_elev,3);
    vyrhs_cvec   = sum((N_row_points_pages.*vny_pages).* weights_pages .* J_det_points_elev,3);

end
