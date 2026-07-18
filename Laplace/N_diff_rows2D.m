% N_diff_rows2D - Compute derivatives of Lagrange shape functions.
%
% FILE: N_diff_rows2D.m
% DESCRIPTION:
% Build 2D arrays of shape function derivatives (dN/dxi, dN/deta) for
% tensor-product Lagrange elements by combining 1D derivatives in xi and eta.
%
% Inputs:
%   K (variable): Number of nodes per element side for the solution element (K = order+1).
%   ind_i (variable): Row index array for tensor-product Lagrange shape functions.
%   ind_j (variable): Column index array for tensor-product Lagrange shape functions.
%   naturalCoords_row_pages (variable): Vector of physical coordinates.
% Outputs:
%   none
function N_diff_rows_pages=N_diff_rows2D(K,ind_i,ind_j,naturalCoords_row_pages)
    xi_elev=naturalCoords_row_pages(1,1,:);
    eta_elev=naturalCoords_row_pages(1,2,:);

    % The general equation is:
    % dN_xi = dN_i(xi).*N_j(eta),
    % dN_eta=N_i(xi).*dN_j(eta)
    Num_nodes_1D=round(sqrt(K));
    N_xi_row_pages=N_row(Num_nodes_1D,xi_elev);
    N_eta_row_pages=N_row(Num_nodes_1D,eta_elev);

    N_xi_diff_row_pages=N_diff_rows(Num_nodes_1D,xi_elev);
    N_eta_diff_row_pages=N_diff_rows(Num_nodes_1D,eta_elev);

    N_diff_rows_pages=[N_xi_diff_row_pages(1,ind_i,:).*N_eta_row_pages(1,ind_j,:);N_xi_row_pages(1,ind_i,:).*N_eta_diff_row_pages(1,ind_j,:)];
end
