% N_diff_rows - Compute derivatives of Lagrange shape functions.
%
% FILE: N_diff_rows.m
% DESCRIPTION:
% Compute 1D derivatives of Lagrange shape functions at given natural
% coordinate points. Used as a building block for 2D tensor-product shapes.
%
% Inputs:
%   K (variable): Number of nodes per element side for the solution element (K = order+1).
%   xi_elev (variable): Natural coordinate evaluation points along the xi direction.
% Outputs:
%   none
function N_diff_row_pages=N_diff_rows(K,xi_elev)
    N_points=size(xi_elev,3); % number of points at which shape functions are to be computed
    xi_i_row=linspace(-1,1,K);
    xi_i_row=xi_i_row([1,K,2:K-1]); % reordering nodes, because node 2 is always at the end
    tmp1=(xi_elev-xi_i_row.')./(xi_i_row-xi_i_row.');
    tmp1((1:K+1:K^2)+(0:N_points-1).'*K^2)=[];
    tmp1=reshape(tmp1,K-1,K,[]);
    if(K==2)
        ind_j_vec=[];
    else
        ind_j_vec=flip(nchoosek(1:K-1,K-2),1).'; % example: for K=4 -> j=[[2,3],[1,3],[1,2]]
    end
    ind_j_vec=reshape(ind_j_vec,[],1,K-1); % reshape into 3rd dimension, because the counter j is the 3rd nested counter in equation 5.68
    tmp1=reshape(tmp1(ind_j_vec,:,:),[],K,N_points,K-1);
    tmp1=prod(tmp1,1);
    tmp1=reshape(tmp1,K-1,K,N_points); % this is a matrix carrying the product terms, i.e. (xi_elev-xi_i_row(j))/(xi_i_row(k)-xi_i_row(j))
    tmp2=1./(xi_i_row-xi_i_row.'); % i.e. 1/(xi_i_row(k)-xi_i_row(i))
    tmp2((1:K+1:K^2))=[];
    tmp2=reshape(tmp2,K-1,K,[]); % this is a matrix carrying the coefficients of the product terms
    N_diff_row_pages=sum(tmp2.*tmp1,1);
end
