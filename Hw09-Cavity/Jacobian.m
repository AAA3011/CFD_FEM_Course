% Jacobian - Compute Jacobian and map natural to physical derivatives.
%
%
% Computes Jacobian matrices and maps shape-function derivatives from
% natural coordinates to physical coordinates at integration points.
%
% Inputs:
%   ElementPhysCoords_rows (variable):Physical element nodal coordinates arranged row-wise.
%   N_diff_NatCoords_points_rows_pages (variable): Vector of physical coordinates.
% Outputs:
%   J_det_points_elev : Jacobian determinant at each Gauss integration point.
%   N_diff_PhysCoords_points_rows_pages : Physical derivatives of shape functions at each Gauss point.
function [J_det_points_elev,N_diff_PhysCoords_points_rows_pages]=Jacobian(ElementPhysCoords_rows,N_diff_NatCoords_points_rows_pages)

% Outputs: J_det_points_elev, N_diff_PhysCoords_points_rows_pages.

    N_pages=size(N_diff_NatCoords_points_rows_pages,3);
    
    J_det_points_elev=nan(1,1,N_pages);
    
    ElementPhysCoords_pages             = repmat(ElementPhysCoords_rows,1,1,N_pages);
    J_points_pages                      = pagemtimes(N_diff_NatCoords_points_rows_pages,ElementPhysCoords_pages);
    N_diff_PhysCoords_points_rows_pages = pagemldivide(J_points_pages, N_diff_NatCoords_points_rows_pages);
    
    for n=1:N_pages
        J_det_points_elev(n)=det(J_points_pages(:,:,n));
    end
end

