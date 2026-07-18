% main - Main simulation driver.
%
% It creates the mesh, builds the global matrices, solves the pressure and
% velocity fields, plots the results, and saves the solution to a .mat file.
%
% Syntax: main()
%
% Inputs: none
%
% Outputs: none
%

clc;clearvars;close all;
tic
%% Global Parameters
epslon    = 0.001;
Re        = 10;
timeStep  = 0.0005;

%% Mesh setup
numElementsX = 4;
numElementsY = 4;
numNodesX    = numElementsX + 1;
numNodesY    = numElementsY + 1;
totNumNodes  = numNodesX * numNodesY;
lengthX      = 1;
lengthY      = 1;
xMin         = 0;
yMin         = 0;

[connectivityMatrix_mat,xCoord_vec,yCoord_vec] = structuredMesh(numElementsX,numElementsY,lengthX,lengthY,xMin,yMin);

%% Precompute geometric data for all elements
numGaussPoints = 4;
n_ElementType  = 1;
numElements    = size(connectivityMatrix_mat,1);
elementData    = cell(numElements, 1);

for elementNumber = 1:numElements

    elementData{elementNumber}.Nodes           = connectivityMatrix_mat(elementNumber,:);
    elementData{elementNumber}.xNodesVals_vec  = xCoord_vec(elementData{elementNumber}.Nodes);
    elementData{elementNumber}.yNodesVals_vec  = yCoord_vec(elementData{elementNumber}.Nodes);

    [weights_pages,N_row_points_pages,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,ElementPhysCoords_mat] = getGaussRelated(elementData{elementNumber}.xNodesVals_vec,elementData{elementNumber}.yNodesVals_vec,n_ElementType,numGaussPoints);

    elementData{elementNumber}.weights    = weights_pages;
    elementData{elementNumber}.N_row      = N_row_points_pages;
    elementData{elementNumber}.J_det      = J_det_points_elev;
    elementData{elementNumber}.dNdX       = N_diff_PhysCoords_points_rows_pages(1, :, :);
    elementData{elementNumber}.dNdY       = N_diff_PhysCoords_points_rows_pages(2, :, :);
    elementData{elementNumber}.physCoords = ElementPhysCoords_mat;

end

%% Build global matrices
pK_mat   = zeros(totNumNodes,totNumNodes);
uMLV_col    = zeros(totNumNodes,1);
vMLV_col    = zeros(totNumNodes,1);

for elementNumber = 1:numElements

    elementNodes_vec = elementData{elementNumber}.Nodes;

    [pk_mat,umlv_col,vmlv_col] = computeLocalNSP(elementData{elementNumber});

    % Generate grid indices matching the local matrix structure
    [R, C] = ndgrid(elementNodes_vec, elementNodes_vec);
    
    % Store into triplet vectors
    idx = cnt + (1:numNodesPerEle^2);
    I(idx) = R(:);
    J(idx) = C(:);
    V(idx) = pk_mat(:);
    cnt    = cnt + numNodesPerEle^2;
    
    % Accumulate column vectors
    uMLV_col(elementNodes_vec)  = uMLV_col(elementNodes_vec) + umlv_col';
    vMLV_col(elementNodes_vec)  = vMLV_col(elementNodes_vec) + vmlv_col';
end

% Create the Global Sparse Matrix directly without high memory overhead
pK_mat = sparse(I, J, V, totNumNodes, totNumNodes);


%% Boundary condition locations
ind_vec = 1:length(xCoord_vec);

pBCCondition_vec       = (xCoord_vec == lengthX/2) & (yCoord_vec == 0);
upperWallNoCorners_vec = (yCoord_vec == lengthY)   & (xCoord_vec ~= lengthX) & (xCoord_vec ~= 0);
rightWallNoCorner_vec  = (xCoord_vec == lengthX)   & (yCoord_vec ~= lengthY);
leftWallNoCorner_vec   = (xCoord_vec == 0)         & (yCoord_vec ~= lengthY);
upperLeftCorner_vec    = (xCoord_vec == 0)         & (yCoord_vec == lengthY);
upperRightCorner_vec   = (xCoord_vec == lengthX)   & (yCoord_vec == lengthY);
lowerWall_vec          =  yCoord_vec == 0;
xEqualHalf             =  xCoord_vec == 0.5;

%% Prepare the pressure matrix
pBCNodes_col                      = ind_vec(pBCCondition_vec);
pK_mat(pBCNodes_col,:)            = 0;
pK_mat(pBCNodes_col,pBCNodes_col) = eye(length(pBCNodes_col));

%% Initial values
un_col = zeros(totNumNodes,1);
vn_col = zeros(totNumNodes,1);

%% Iteration controls
itr     = 0;
totitrs = 1130;
pK_mat  = decomposition(pK_mat);   % this line to enhance performance

while(itr < totitrs)
    pressureSolution_col = computePressure(elementData,numElements,pBCNodes_col,pK_mat,totNumNodes,un_col,vn_col,epslon);

    uRHS_col    = zeros(totNumNodes,1);
    vRHS_col    = zeros(totNumNodes,1);
    ReLocal_col = zeros(numElements,1);
    CFL_col     = zeros(numElements,1);
    xCenter_col = zeros(numElements,1);
    yCenter_col = zeros(numElements,1);

    for elementNumber = 1:numElements

        elementNodes_vec           = elementData{elementNumber}.Nodes;

        unNodes_col                = un_col(elementNodes_vec,1);
        vnNodes_col                = vn_col(elementNodes_vec,1);
        pressureSolutionNodes_col  = pressureSolution_col(elementNodes_vec,1);

        [urhs_col,vrhs_col,ReLocal,cfl,xCenter,yCenter] = computeLocalNS(elementData{elementNumber},unNodes_col,vnNodes_col,pressureSolutionNodes_col,Re,timeStep);

        uRHS_col(elementNodes_vec)  = uRHS_col(elementNodes_vec)  + urhs_col;
        vRHS_col(elementNodes_vec)  = vRHS_col(elementNodes_vec)  + vrhs_col;
        ReLocal_col(elementNumber)  = ReLocal;
        CFL_col(elementNumber)      = cfl;
        xCenter_col(elementNumber)  = xCenter;
        yCenter_col(elementNumber)  = yCenter;

    end

    %% Apply boundary conditions for u velocity
    uRHS_col(ind_vec(lowerWall_vec))          = 0;
    uRHS_col(ind_vec(upperWallNoCorners_vec)) = 1;
    uRHS_col(ind_vec(rightWallNoCorner_vec))  = 0;
    uRHS_col(ind_vec(leftWallNoCorner_vec))   = 0;
    uRHS_col(ind_vec(upperLeftCorner_vec))    = 1;
    uRHS_col(ind_vec(upperRightCorner_vec))   = 1;

    uMLV_col(ind_vec(lowerWall_vec))          = 1;
    uMLV_col(ind_vec(upperWallNoCorners_vec)) = 1;
    uMLV_col(ind_vec(rightWallNoCorner_vec))  = 1;
    uMLV_col(ind_vec(leftWallNoCorner_vec))   = 1;
    uMLV_col(ind_vec(upperLeftCorner_vec))    = 1;
    uMLV_col(ind_vec(upperRightCorner_vec))   = 1;

    %% Apply boundary conditions for v velocity
    vRHS_col(ind_vec(lowerWall_vec))          = 0;
    vRHS_col(ind_vec(upperWallNoCorners_vec)) = 0;
    vRHS_col(ind_vec(rightWallNoCorner_vec))  = 0;
    vRHS_col(ind_vec(leftWallNoCorner_vec))   = 0;
    vRHS_col(ind_vec(upperLeftCorner_vec))    = 0;
    vRHS_col(ind_vec(upperRightCorner_vec))   = 0;

    vMLV_col(ind_vec(lowerWall_vec))          = 1;
    vMLV_col(ind_vec(upperWallNoCorners_vec)) = 1;
    vMLV_col(ind_vec(rightWallNoCorner_vec))  = 1;
    vMLV_col(ind_vec(leftWallNoCorner_vec))   = 1;
    vMLV_col(ind_vec(upperLeftCorner_vec))    = 1;
    vMLV_col(ind_vec(upperRightCorner_vec))   = 1;

    %% Solve for u velocity
    uSolution_col = uRHS_col./uMLV_col;

    %% Solve for v velocity
    vSolution_col = vRHS_col./vMLV_col;

    %% Use the new values as the next initial guess
    un_col  = uSolution_col;
    vn_col  = vSolution_col;

    %% Step forward in time
    itr = itr + 1;
    disp(itr)
end
toc

%% Save results
% save('cavity_resultsRe10.mat', 'uSolution_col', 'vSolution_col', 'pressureSolution_col', 'xCoord_vec', 'yCoord_vec', 'connectivityMatrix_mat');
