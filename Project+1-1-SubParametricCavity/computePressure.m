% computePressure - Compute the global pressure system.
%
% Description: Build pressure RHS from velocity divergence and solve for
% nodal pressure values using preassembled pressure matrix.
% Inputs:
%   elementData (variable): Struct with precomputed element integration data (weights, shape functions, derivatives, physCoords)
%   connectivityMatrix_mat (matrix): Element connectivity matrix (nElements x nodesPerElement)
%   pBCNodes_col (column vector): Indices of nodes with pressure Dirichlet BCs
%   pK_mat (matrix): Global pressure matrix (sparse)
%   totNumNodes (variable): Total number of nodes (scalar)
%   un_col (column vector): Nodal u-velocity column vector (N x 1)
%   vn_col (column vector): Nodal v-velocity column vector (N x 1)
%   epslon (variable): Small penalty/stabilization parameter used in the pressure Poisson RHS.
% Outputs:
%   pressureSolution_col : Nodal pressure column vector (N x 1)
function[pressureSolution_col] = computePressure(elementData,connectivityMatrix_mat,pBCNodes_col,pK_mat,totNumNodes,un_col,vn_col,epslon)

%% pRHS initialization
pRHS_col = zeros(totNumNodes,1);

for eleNum = 1:size(connectivityMatrix_mat,1)

    elementNodes_vec   = elementData{eleNum}.Nodes;
    weights_pages      = elementData{eleNum}.weights;
    N_row_points_pages = elementData{eleNum}.N_row;
    J_det_points_elev  = elementData{eleNum}.J_det;
    dNdX_vec_pages     = elementData{eleNum}.dNdX;
    dNdY_vec_pages     = elementData{eleNum}.dNdY;

    unNodes_col = un_col(elementNodes_vec,1);
    vnNodes_col = vn_col(elementNodes_vec,1);

    unx_pages   = pagemtimes(dNdX_vec_pages,unNodes_col);
    vny_pages   = pagemtimes(dNdY_vec_pages,vnNodes_col);

    prhs_col  = sum(((-1/epslon).*(unx_pages+vny_pages).*permute(N_row_points_pages,[2 1 3]))  .* weights_pages .* J_det_points_elev,3);

    pRHS_col(elementNodes_vec) = pRHS_col(elementNodes_vec) + prhs_col;
end

%% Pressure Boundary Conditions
pRHS_col(pBCNodes_col)  = 0;

%% Pressure Solution
pressureSolution_col = pK_mat \ pRHS_col;
end
