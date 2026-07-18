% getGaussRelated - Return Gauss Related.
%
% FILE: getGaussRelated.m
% DESCRIPTION:
% Computes shape functions, derivatives (natural and physical coords),
% Jacobian determinants, and quadrature weights at all Gauss points for
% a single element. Central hub for precomputing element geometry data.
%
% Inputs:
%   xNodesValsGeo_vec (column vector): Element nodal x-coordinates.
%   yNodesValsGeo_vec (column vector): Element nodal y-coordinates.
%   K (variable): Number of nodes in one element direction for the current element type.
%   numGaussPoints (variable): Number of Gauss points used per dimension
% Outputs:
%   weights_pages : Gauss quadrature weights at each integration point, stored page-wise.
%   N_row_points_pages : Shape function values evaluated at Gauss points, stored page-wise.
%   J_det_points_elev : Jacobian determinant at each Gauss integration point.
%   N_diff_PhysCoords_points_rows_pages : Physical derivatives of shape functions at each Gauss point.
%   ElementGeoCoords_mat : Physical element nodal coordinates matrix.
function [weights_pages,N_row_points_pages,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,ElementGeoCoords_mat] = getGaussRelated(xNodesValsGeo_vec,yNodesValsGeo_vec,K,numGaussPoints)

    [gaussPoints_col, gaussWeights_col] = gaussLegendreQuad(numGaussPoints);

    [xi_mat, eta_mat] = meshgrid(gaussPoints_col, gaussPoints_col);    % make all different combinations of xi , eta
    coords_row_pages  = reshape([xi_mat(:), eta_mat(:)]', 1, 2, []);   % So now we have the coordinates for all compinations stored page wise
    weights_mat       = gaussWeights_col(:) * gaussWeights_col(:)';    % Make a weights tensor which means multibled weights for each coordinates set
    weights_row       = weights_mat(:)';                               % So now we have the cooresponding weights stored in a row
    weights_pages     = reshape(weights_row,1,1,[]);

    % Compute natural shape functions and their derivatives
    ElementGeoCoords_mat = [xNodesValsGeo_vec; yNodesValsGeo_vec]';
    [N_row_points_pages, N_diff_rows_points_pages] = shapeFunctions(K, coords_row_pages);

    % Compute Jacobian and physical derivatives
    [J_det_points_elev, N_diff_PhysCoords_points_rows_pages] = Jacobian(ElementGeoCoords_mat, N_diff_rows_points_pages);
end
