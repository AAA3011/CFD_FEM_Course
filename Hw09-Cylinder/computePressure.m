% computePressure - Compute the global pressure system.
%
% FILE: computePressure.m
% DESCRIPTION:
% Solves pressure Poisson equation at each time step. Assembles RHS from
% velocity divergence, applies pressure boundary conditions, and solves
% the linear system for nodal pressure values.
%
% Inputs:
%   elementData (variable): Struct with precomputed element integration data (weights, shape functions, derivatives, physCoords)
%   numElements (variable): Number of elements (scalar)
%   pBCNodes_vec (column vector): Indices of nodes with pressure Dirichlet BCs
%   pK_mat (matrix): Global pressure matrix (sparse)
%   totNumNodes (variable): Total number of nodes (scalar)
%   un_col (column vector): Nodal u-velocity column vector (N x 1)
%   vn_col (column vector): Nodal v-velocity column vector (N x 1)
%   epslon (variable): Small penalty/stabilization parameter used in the pressure Poisson RHS.
% Outputs:
%   pressureSolution_col : Nodal pressure column vector (N x 1)
function[pressureSolution_col] = computePressure(elementData,numElements,pBCNodes_vec,pK_mat,totNumNodes,un_col,vn_col,epslon)

    %% Initialize the pressure source term
    pRHS_col = zeros(totNumNodes,1);

    %% Loop over elements
    for eleNum = 1:numElements
    
        elementNodes_vec = elementData{eleNum}.Nodes;
        unNodes_col     = un_col(elementNodes_vec,1);
        vnNodes_col     = vn_col(elementNodes_vec,1);
    

        dNdX_vec_pages     = elementData{eleNum}.dNdX;
        dNdY_vec_pages     = elementData{eleNum}.dNdY;
        weights_pages      = elementData{eleNum}.weights;
        J_det_points_elev  = elementData{eleNum}.J_det;
        N_row_points_pages = elementData{eleNum}.N_row;


        unx_pages = pagemtimes(dNdX_vec_pages,unNodes_col);
        vny_pages = pagemtimes(dNdY_vec_pages,vnNodes_col);

        
        prhs_col  = sum(((-1/epslon).*(unx_pages+vny_pages).*permute(N_row_points_pages,[2 1 3]))  .* weights_pages .* J_det_points_elev,3);
        
        pRHS_col(elementNodes_vec) = pRHS_col(elementNodes_vec) + prhs_col;
    end
    
    %% Apply pressure boundary conditions
    pRHS_col(pBCNodes_vec)  = 0;
    
    %% Solve the pressure system
    pressureSolution_col = pK_mat \ pRHS_col;
end
