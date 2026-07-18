function [Diver,Ku_u_mat_pages,Ku_v_mat_pages,Ku_p_mat_pages,Kv_u_mat_pages,Kv_v_mat_pages,Kv_p_mat_pages,Kp_u_mat_pages,Kp_v_mat_pages,Kp_p_mat_pages,Ru_col_pages,Rv_col_pages,Rp_col_pages] = computeLocal(geomData,ukNodes_cvec,vkNodes_cvec,pkNodes_cvec,timeStep,unNodes_cvec,vnNodes_cvec,Re,B)

    %% Extract Element Informations
    weights_pages         = geomData.weights;
    N_row_points_pages    = geomData.N_row;
    J_det_points_elev     = geomData.J_det;
    dNdX_vec_pages        = geomData.dNdX;
    dNdY_vec_pages        = geomData.dNdY;
    ElementPhysCoords_mat = geomData.physCoords;
    
    %% Coordinates
    xphysCoord_vec = ElementPhysCoords_mat(1,:);
    yphysCoord_vec = ElementPhysCoords_mat(2,:);
    
    %% interpolation For Values per Element inprevious step
    ukx_pages  = pagemtimes(dNdX_vec_pages,ukNodes_cvec);
    uky_pages  = pagemtimes(dNdY_vec_pages,ukNodes_cvec);
    vky_pages  = pagemtimes(dNdY_vec_pages,vkNodes_cvec);
    vkx_pages  = pagemtimes(dNdX_vec_pages,vkNodes_cvec);
    uk_pages   = pagemtimes(N_row_points_pages,ukNodes_cvec);
    un_pages   = pagemtimes(N_row_points_pages,unNodes_cvec);
    vk_pages   = pagemtimes(N_row_points_pages,vkNodes_cvec);
    vn_pages   = pagemtimes(N_row_points_pages,vnNodes_cvec);
    pkx_pages  = pagemtimes(dNdX_vec_pages,pkNodes_cvec);
    pky_pages  = pagemtimes(dNdY_vec_pages,pkNodes_cvec);
    unx_pages  = pagemtimes(dNdX_vec_pages,unNodes_cvec);
    vny_pages  = pagemtimes(dNdY_vec_pages,vnNodes_cvec);
    %% Parameteres and Weight calculations
    le1 = sqrt((yphysCoord_vec(3)-yphysCoord_vec(1))^2+(xphysCoord_vec(3)-xphysCoord_vec(1))^2);
    le2 = sqrt((yphysCoord_vec(4)-yphysCoord_vec(2))^2+(xphysCoord_vec(4)-xphysCoord_vec(2))^2);
    elementLength  = max(le1,le2);
    
    Ve_pages      = sqrt(uk_pages.^2 + vk_pages.^2);
    thu           = 1./sqrt((2/(timeStep))^2 + ((2*Ve_pages)/elementLength).^2 + (4/(elementLength^2*Re)).^2);
    W_cvec_pages  = permute(N_row_points_pages,[2 1 3]) + thu.*(uk_pages.*(permute(dNdX_vec_pages,[2 1 3])) + vk_pages.*(permute(dNdY_vec_pages,[2 1 3])));
    
    %% Preparing NNx, NNy, NxNx and NyNy
    NiNj_mat_pages   = pagemtimes(permute(N_row_points_pages,[2 1 3]), N_row_points_pages);
    NiNjx_mat_pages  = pagemtimes(permute(N_row_points_pages,[2 1 3]), dNdX_vec_pages);
    NiNjy_mat_pages  = pagemtimes(permute(N_row_points_pages,[2 1 3]), dNdY_vec_pages);
    NixNjx_mat_pages = pagemtimes(permute(dNdX_vec_pages,[2 1 3]),dNdX_vec_pages);
    NiyNjy_mat_pages = pagemtimes(permute(dNdY_vec_pages,[2 1 3]),dNdY_vec_pages);
    
    %% Preparing NxW, NW, NyW
    NixWj  = pagemtimes(W_cvec_pages, dNdX_vec_pages);
    NiWj   = pagemtimes(W_cvec_pages, N_row_points_pages);
    NiyWj  = pagemtimes(W_cvec_pages, dNdY_vec_pages);
    
    %% Pressure-Poisson Equation Contribution
    Kp_u_mat_pages = -2*sum((ukx_pages .* NiNjx_mat_pages + vkx_pages .* NiNjy_mat_pages) .* weights_pages .* J_det_points_elev, 3);
    Kp_v_mat_pages = -2*sum((vky_pages .* NiNjy_mat_pages + uky_pages .* NiNjx_mat_pages) .* weights_pages .* J_det_points_elev, 3);
    Kp_p_mat_pages = sum((NixNjx_mat_pages + NiyNjy_mat_pages) .* weights_pages .* J_det_points_elev, 3);

    Diver_be = (unx_pages + vny_pages);
    % termRhs1Rp   = (ukx_pages.^2 + 2*uky_pages.*vkx_pages + vky_pages.^2 - (1/timeStep)*Diver_be) .* permute(N_row_points_pages,[2 1 3]);
    termRhs1Rp   = (ukx_pages.^2 + 2*uky_pages.*vkx_pages + vky_pages.^2 ) .* permute(N_row_points_pages,[2 1 3]);
    % disp( sum((1/timeStep)*(unx_pages + vny_pages) .* weights_pages .* J_det_points_elev, 3))
    termRhs2Rp   = pkx_pages .* permute(dNdX_vec_pages,[2 1 3]) + pky_pages .* permute(dNdY_vec_pages,[2 1 3]);
    Rp_col_pages = sum((termRhs1Rp - termRhs2Rp) .* weights_pages .* J_det_points_elev, 3);
    Diver = sum(Diver_be .* weights_pages .* J_det_points_elev, 3);
    %% X-Momentum Equation Contribution
    term1kuu = NiNj_mat_pages;
    term2kuu = uk_pages  .* NixWj ;
    term3kuu = ukx_pages .* NiWj  ;
    term4kuu = vk_pages  .* NiyWj ;
    term5kuu = NixNjx_mat_pages + NiyNjy_mat_pages;
    
    term1Ru  = (uk_pages .* ukx_pages + vk_pages .* uky_pages).* W_cvec_pages;
    term2Ru  = ukx_pages .* permute(dNdX_vec_pages,[2 1 3]) + uky_pages .* permute(dNdY_vec_pages,[2 1 3]);
    term3Ru  = (pkx_pages.* W_cvec_pages) + ((uk_pages - un_pages)/(timeStep)).*permute(N_row_points_pages,[2 1 3]);
    
    Ku_u_mat_pages = sum((((1/timeStep) * term1kuu) + term2kuu + term3kuu + term4kuu + (B*term5kuu)).* weights_pages .* J_det_points_elev,3);
    Ku_v_mat_pages = sum((uky_pages .* NiWj ).* weights_pages .* J_det_points_elev,3);
    Ku_p_mat_pages = sum((NixWj).* weights_pages .* J_det_points_elev,3);
    
    Ru_col_pages = sum(( -term1Ru  - B * term2Ru - term3Ru).* weights_pages .* J_det_points_elev,3);
    
    %% Y-Momentum Equation Contribution
    term1kvv = NiNj_mat_pages;
    term2kvv = uk_pages  .* NixWj;
    term3kvv = vky_pages .* NiWj;
    term4kvv = vk_pages  .* NiyWj;
    term5kvv = NixNjx_mat_pages + NiyNjy_mat_pages;
    
    term1Rv = (uk_pages  .* vkx_pages + vk_pages .* vky_pages) .* W_cvec_pages;
    term2Rv = (vkx_pages .* permute(dNdX_vec_pages,[2 1 3]) + vky_pages .* permute(dNdY_vec_pages,[2 1 3]));
    term3Rv = pky_pages.* W_cvec_pages + ((vk_pages - vn_pages)/timeStep) .* permute(N_row_points_pages,[2 1 3]);
    
    Kv_u_mat_pages = sum((vkx_pages .* NiWj).* weights_pages .* J_det_points_elev,3);
    Kv_v_mat_pages = sum(((1/timeStep).* term1kvv + term2kvv + term3kvv + term4kvv + B.* term5kvv).* weights_pages .* J_det_points_elev,3);
    Kv_p_mat_pages = sum((NiyWj).* weights_pages .* J_det_points_elev,3);
    
    Rv_col_pages = sum((-term1Rv - B*term2Rv - term3Rv).* weights_pages .* J_det_points_elev,3);

end
