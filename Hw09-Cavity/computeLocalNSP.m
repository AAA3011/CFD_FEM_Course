% computeLocalNSP - Compute local NSP.
%
%
% Assembles the element pressure stiffness (Poisson) contribution and
% returns auxiliary mass-like vectors used during global assembly.
%
% Inputs:
%   elementData (variable): Struct with precomputed element integration data (weights, shape functions, derivatives, physCoords)
% Outputs:
%   pk_mat : Global pressure matrix (sparse)
%   umlv_col : Local u-momentum mass-like contribution vector.
%   vmlv_col : Local v-momentum mass-like contribution vector.
function [pk_mat,umlv_col,vmlv_col] = computeLocalNSP(elementData)
% Outputs: pk_mat, umlv_col, vmlv_col.


    weights_pages         = elementData.weights;
    N_row_points_pages    = elementData.N_row;
    J_det_points_elev     = elementData.J_det;
    dNdX_vec_pages        = elementData.dNdX;
    dNdY_vec_pages        = elementData.dNdY;

    pk_mat      = sum((pagemtimes(permute(dNdX_vec_pages,[2 1 3]),dNdX_vec_pages) + pagemtimes(permute(dNdY_vec_pages,[2 1 3]),dNdY_vec_pages))  .* weights_pages .* J_det_points_elev,3);
    umlv_col   = sum(N_row_points_pages .* weights_pages .* J_det_points_elev,3);
    vmlv_col   = sum(N_row_points_pages .* weights_pages .* J_det_points_elev,3);
end
