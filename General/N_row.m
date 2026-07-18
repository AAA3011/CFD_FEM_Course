% N_row - Compute 1D/2D Lagrange shape function values.
%
% FILE: N_row.m
% DESCRIPTION:
% Compute 1D Lagrange shape function values at given natural coordinate
% points for a set of nodes placed on the reference interval [-1,1].
%
% Inputs:
%   K (variable): Number of nodes per element side for the solution element (K = order+1).
%   xi_elev (variable): Natural coordinate evaluation points along the xi direction.
% Outputs:
%   none
function N_row_pages=N_row(K,xi_elev)

    xi_i_row=linspace(-1,1,K);
    xi_i_row=xi_i_row([1,K,2:K-1]); % reordering nodes, because node 2 is always at the end
    N_row_tmp=(xi_elev-xi_i_row.')./(xi_i_row-xi_i_row.');
    ind_vec=(1:K+1:K^2)+(0:size(xi_elev,3)-1).'*K^2;
    N_row_tmp(ind_vec)=[];
    N_row_pages=prod(reshape(N_row_tmp,K-1,K,[]),1);
end
