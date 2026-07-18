% Jacobian - Compute Jacobian and map natural to physical derivatives.

% Jacobian - Compute Jacobian and map natural to physical derivatives.
%
% Description: Compute Jacobian matrices/determinants and transform
% shape derivatives from natural to physical coordinates.
% Inputs:
%   ElementPhysCoords_rows (variable):Physical element nodal coordinates arranged row-wise.
%   N_diff_NatCoords_points_rows_pages (variable): Vector of physical coordinates.
%   N_diff_Geo_rows_points_pages (variable): Geometry shape function derivatives at geometry Gauss points.
% Outputs:
%   J_det_points_elev : Jacobian determinant at each Gauss integration point.
%   N_diff_PhysCoords_points_rows_pages : Physical derivatives of shape functions at each Gauss point.
function [J_det_points_elev,N_diff_PhysCoords_points_rows_pages]=Jacobian(ElementPhysCoords_rows,N_diff_NatCoords_points_rows_pages,N_diff_Geo_rows_points_pages)

    N_pages    = size(N_diff_NatCoords_points_rows_pages,3);
    NGeo_pages = size(N_diff_Geo_rows_points_pages,3);
    
    J_det_points_elev = zeros(1,1,N_pages);

    N_diff_Geo_rows_points_pages         = permute(N_diff_Geo_rows_points_pages,[2 1 3]);
    ElementGeoCoords_pages               = repmat(ElementPhysCoords_rows,1,1,NGeo_pages);
    J_points_pages                       = pagetranspose(pagemtimes(ElementGeoCoords_pages,N_diff_Geo_rows_points_pages));
    N_diff_PhysCoords_points_rows_pages  = pagemldivide(J_points_pages, N_diff_NatCoords_points_rows_pages);
    
    for n=1:N_pages
        J_det_points_elev(1,1,n) = det(J_points_pages(:,:,n));
    end
end

