% computeLocalNS - Compute local NS.
%
% Description: Compute local Navier-Stokes residual contributions and
% auxiliary quantities for one element using Gauss quadrature.
% Inputs:
%   elementData (variable): Struct with precomputed element integration data (weights, shape functions, derivatives, physCoords)
%   unNodes_col : Nodal u-velocity values for the current element.
%   vnNodes_col : Nodal v-velocity values for the current element.
%   pressureSolutionNodes_col : Nodal pressure values for the current element.
%   Re (variable): Reynolds number (scalar)
%   timeStep (variable): Time step size (scalar)
% Outputs:
%   urhs_col : Local element u-momentum RHS vector (N x 1).
%   vrhs_col : Local element v-momentum RHS vector (N x 1).
%   ReLocal : Reynolds number (scalar)
%   cfl : Local CFL number for the element (dimensionless).
%   xCenter : Element centroid x-coordinate.
%   yCenter : Element centroid y-coordinate.
%   uxrhs_col : Nodal u-velocity values.
%   uyrhs_col : Nodal u-velocity values.
%   vxrhs_col : Nodal v-velocity values.
%   vyrhs_col : Nodal v-velocity values.
function [urhs_col,vrhs_col,ReLocal,cfl,xCenter,yCenter,uxrhs_col,uyrhs_col,vxrhs_col,vyrhs_col] = computeLocalNS(elementData,unNodes_col,vnNodes_col,pressureSolutionNodes_col,Re,timeStep)
 
    xphysCoord_vec     = elementData.xNodesGVals_vec;
    yphysCoord_vec     = elementData.yNodesGVals_vec;

    weights_pages      = elementData.weights;
    N_row_points_pages = elementData.N_row;
    J_det_points_elev  = elementData.J_det;
    dNdX_vec_pages     = elementData.dNdX;
    dNdY_vec_pages     = elementData.dNdY;

    unx_pages  = pagemtimes(dNdX_vec_pages,unNodes_col);
    uny_pages  = pagemtimes(dNdY_vec_pages,unNodes_col);
    vny_pages  = pagemtimes(dNdY_vec_pages,vnNodes_col);
    vnx_pages  = pagemtimes(dNdX_vec_pages,vnNodes_col);
    un_pages   = pagemtimes(N_row_points_pages,unNodes_col);
    vn_pages   = pagemtimes(N_row_points_pages,vnNodes_col);
    pnx_pages  = pagemtimes(dNdX_vec_pages,pressureSolutionNodes_col);
    pny_pages  = pagemtimes(dNdY_vec_pages,pressureSolutionNodes_col);

    le1            = sqrt((yphysCoord_vec(3)-yphysCoord_vec(1))^2+(xphysCoord_vec(3)-xphysCoord_vec(1))^2);
    le2            = sqrt((yphysCoord_vec(4)-yphysCoord_vec(2))^2+(xphysCoord_vec(4)-xphysCoord_vec(2))^2);
    elementLength  = max(le1,le2);

    Ve_pages    = sqrt(un_pages.^2 + vn_pages.^2);
    thu         = 1./sqrt((2/(timeStep))^2 + ((2*Ve_pages)/elementLength).^2 + (4./((elementLength^2)*Re)^2));
    W_col_pages = permute(N_row_points_pages,[2 1 3]) + thu.*(un_pages.*(permute(dNdX_vec_pages,[2 1 3])) + vn_pages.*(permute(dNdY_vec_pages,[2 1 3])));

    u           = sum(unNodes_col)/4;
    v           = sum(vnNodes_col)/4;

    uCenter     = sqrt((u)^2 + (v)^2);
    ReLocal     = uCenter * Re * elementLength;
    cfl         = uCenter * timeStep / elementLength;
    xCenter     = sum(xphysCoord_vec)/4;
    yCenter     = sum(yphysCoord_vec)/4;
    urhs_col    = sum(((unNodes_col.*permute(N_row_points_pages,[2 1 3]))- (timeStep * (un_pages.*unx_pages + vn_pages.*uny_pages) .* W_col_pages) - (timeStep * pnx_pages .* permute(N_row_points_pages,[2 1 3])) - ((timeStep/Re)*(unx_pages.*(permute(dNdX_vec_pages,[2 1 3])) + uny_pages.*(permute(dNdY_vec_pages,[2 1 3]))))) .* weights_pages .* J_det_points_elev,3);
    vrhs_col    = sum(((vnNodes_col.*permute(N_row_points_pages,[2 1 3]))- (timeStep * (un_pages.*vnx_pages + vn_pages.*vny_pages) .* W_col_pages) - (timeStep * pny_pages .* permute(N_row_points_pages,[2 1 3])) - ((timeStep/Re)*(vnx_pages.*(permute(dNdX_vec_pages,[2 1 3])) + vny_pages.*(permute(dNdY_vec_pages,[2 1 3]))))) .* weights_pages .* J_det_points_elev,3);
    uxrhs_col   = sum((N_row_points_pages.*unx_pages).* weights_pages .* J_det_points_elev,3);
    uyrhs_col   = sum((N_row_points_pages.*uny_pages).* weights_pages .* J_det_points_elev,3);
    vxrhs_col   = sum((N_row_points_pages.*vnx_pages).* weights_pages .* J_det_points_elev,3);
    vyrhs_col   = sum((N_row_points_pages.*vny_pages).* weights_pages .* J_det_points_elev,3);
    % Div       = sum((unx_pages + vny_pages).* weights_pages .* J_det_points_elev,3);
end
