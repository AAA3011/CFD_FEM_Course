function [J_det_points_elev,N_diff_PhysCoords_points_rows_pages]=Jacobian(ElementPhysCoords_rows,N_diff_NatCoords_points_rows_pages)

    N_pages=size(N_diff_NatCoords_points_rows_pages,3);
    
    J_det_points_elev=nan(1,1,N_pages);
    
    ElementPhysCoords_pages             = repmat(ElementPhysCoords_rows,1,1,N_pages);
    J_points_pages                      = pagemtimes(N_diff_NatCoords_points_rows_pages,ElementPhysCoords_pages);
    N_diff_PhysCoords_points_rows_pages = pagemldivide(J_points_pages, N_diff_NatCoords_points_rows_pages);
    
    for n=1:N_pages
        J_det_points_elev(n)=det(J_points_pages(:,:,n));
    end
    
    % N_pages=size(N_diff_NatCoords_points_rows_pages,3);
%
% J_det_points_elev=nan(1,1,N_pages);
% N_diff_PhysCoords_points_rows_pages=nan(size(N_diff_NatCoords_points_rows_pages));
% for n=1:N_pages
%     J_points_mat=N_diff_NatCoords_points_rows_pages(:,:,n)*ElementPhysCoords_rows;
%     J_det_points_elev(n)=det(J_points_mat);
%     N_diff_PhysCoords_points_rows_pages(:,:,n)=J_points_mat\N_diff_NatCoords_points_rows_pages(:,:,n);
% end

end

