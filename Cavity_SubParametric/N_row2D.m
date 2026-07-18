function N_row_pages=N_row2D(K,ind_i,ind_j,naturalCoords_row_pages)
    xi_elev=naturalCoords_row_pages(1,1,:);
    eta_elev=naturalCoords_row_pages(1,2,:);

    % The general equation is: N = N_i(xi).*N_j(eta)
    Num_nodes_1D=round(sqrt(K));
    N_xi_row_pages=N_row(Num_nodes_1D,xi_elev);
    N_eta_row_pages=N_row(Num_nodes_1D,eta_elev);
    N_row_pages=N_xi_row_pages(1,ind_i,:).*N_eta_row_pages(1,ind_j,:);

    % temp_mat_pages=pagemtimes(N_xi_row_pages,'transpose',N_eta_row_pages,'none');
    % N_row_pages=temp_mat_pages(ind_i,ind_j,:);
end
