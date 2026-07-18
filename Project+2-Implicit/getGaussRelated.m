% getGaussRelated - Return Gauss Related.
%
% FILE: getGaussRelated.m
% DESCRIPTION:
% Prepare Gauss quadrature-related quantities: weights, shape function
% evaluations at Gauss points, Jacobian determinants and physical
% derivatives needed for element assembly.
%
% Inputs:
%   xNodesVals_vec (column vector): Element nodal x-coordinates.
%   yNodesVals_vec (column vector): Element nodal y-coordinates.
%   K (variable): Number of nodes per element side for the solution element (K = order+1).
%   numGaussPoints (variable): Number of Gauss points used per dimension
%   K_G (variable): Number of nodes per element side for the geometric element mesh (geometry order + 1).
%   numGeometryGaussPoints (variable): Number of Gauss points used for the geometric field.
%   order (variable): Polynomial order of the solution element basis.
%   order_G (variable): Polynomial order of the geometric element basis.
% Outputs:
%   weights_pages : Gauss quadrature weights at each integration point, stored page-wise.
%   N_row_points_pages : Shape function values evaluated at Gauss points, stored page-wise.
%   J_det_points_elev : Jacobian determinant at each Gauss integration point.
%   N_diff_PhysCoords_points_rows_pages : Physical derivatives of shape functions at each Gauss point.
%   ElementPhysCoords_mat : Physical element nodal coordinates matrix.
function [weights_pages,N_row_points_pages,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,ElementPhysCoords_mat] = getGaussRelated(xNodesVals_vec,yNodesVals_vec,K,numGaussPoints,K_G,numGeometryGaussPoints,order,order_G)
    
    %% Solution Field Part
    [gaussPoints_col, gaussWeights_col] = gaussLegendreQuad(numGaussPoints);
    
    [xi_mat, eta_mat] = meshgrid(gaussPoints_col, gaussPoints_col);    % make all different combinations of xi , eta
    coords_row_pages  = reshape([xi_mat(:), eta_mat(:)]', 1, 2, []);   % So now we have the coordinates for all compinations stored page wise
    weights_mat       = gaussWeights_col(:) * gaussWeights_col(:)';    % Make a weights tensor which means multibled weights for each coordinates set
    weights_row       = weights_mat(:)';                               % So now we have the cooresponding weights stored in a row
    weights_pages     = reshape(weights_row,1,1,[]);

    xNodesVals_vec = xNodesVals_vec(:);
    yNodesVals_vec = yNodesVals_vec(:);
    ElementPhysCoords_mat = [xNodesVals_vec'; yNodesVals_vec'];
    
    %% Geometry Field Part
    [gaussPointsGeo_col, ~] = gaussLegendreQuad(numGeometryGaussPoints);
    
    [xiGeo_mat, etaGeo_mat] = meshgrid(gaussPointsGeo_col, gaussPointsGeo_col);
    coordsGeo_row_pages     = reshape([xiGeo_mat(:), etaGeo_mat(:)]', 1, 2, []);
    % weightsGeo_mat          = gaussWeightsGeo_col(:) * gaussWeightsGeo_col(:)';
    % weightsGeo_row          = weightsGeo_mat(:)';
    % weightsGeo_pages        = reshape(weightsGeo_row,1,1,[]);
    
    %% Solution field shape function
    [ind_i, ind_j, ~] = getLagrangeIndices(order);

    N_row_points_pages       = N_row2D(K,ind_i,ind_j,coords_row_pages);
    N_diff_rows_points_pages = N_diff_rows2D(K,ind_i,ind_j,coords_row_pages);
    
    %% Geometry shape function
    [ind_i_G, ind_j_G, ~] = getLagrangeIndices(order_G);
    
    N_diff_Geo_rows_points_pages  = N_diff_rows2D(K_G,ind_i_G,ind_j_G,coordsGeo_row_pages);
    
    %% Jacobian calculation and physical derivatives
    [J_det_points_elev, N_diff_PhysCoords_points_rows_pages] = Jacobian(ElementPhysCoords_mat, N_diff_rows_points_pages,N_diff_Geo_rows_points_pages);
    % area     = sum(J_det_points_elev(:) .* weights_pages(:));
    % area_geo = sum(J_det_points_elev(:) .* weightsGeo_pages(:));

end
