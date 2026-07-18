% main - Main simulation driver.
%
%
% Driver script: sets solver parameters, builds the mesh and element data,
% assembles global matrices, performs time-stepping, writes VTK output,
% and saves results to a .mat file.
%
% Inputs: none (script). Edit parameters inside the file (Re, mesh, order, etc.).
% Outputs: VTK files and a .mat file with `uSolution_col`, `vSolution_col`,
% `pressureSolution_col`, coordinates and connectivity arrays.
%
% Syntax: main()
%
% Inputs: none
%
% Outputs: none
%

clc;clearvars;close all;
tic
%% Element informations
clc;clearvars;close all;
tic
%% Element informations
numGaussPoints         = 10;
numGeometryGaussPoints = numGaussPoints;

D = 2;      % Diemention (here is 2D)

order   = 3;                % Order of Solution Element
K       = (order + 1)^D;    % Solution Element number of nodes

order_G = 1;                % Order of Geometric Element
K_G     = (order_G + 1)^D;  % Geometric Element number of nodes

%% Global Parameters
Re        = 10;
epslon    = 0.001;
timeStep  = 0.0005;
totitrs   = 1130;

vtkOutputInterval = 10;
outputFolder = sprintf('CavityResutsRe%dHigher', Re);
makeVTKdir(outputFolder)

%% Connectivity Matrix and Grid Generation
numElementsX = 10;
numElementsY = 10;
numNodesX    = numElementsX + 1;
numNodesY    = numElementsY + 1;

lengthX      = 1;
lengthY      = 1;
xMin         = 0;
yMin         = 0;

[connectivityMatrixGeo_mat,xCoordGeo_vec,yCoordGeo_vec] = structuredMesh(numElementsX,numElementsY,lengthX,lengthY,xMin,yMin);

%% Converting to higher order mesh 
if order > 1
    [quad_high_conn, xCoord_high, yCoord_high] = quad4_to_quadHigh(connectivityMatrixGeo_mat, xCoordGeo_vec, yCoordGeo_vec, order);
    connectivityMatrix_mat = quad_high_conn;
    xCoord_vec  = xCoord_high;
    yCoord_vec  = yCoord_high;
else
    connectivityMatrix_mat = connectivityMatrixGeo_mat;
    xCoord_vec  = xCoordGeo_vec;
    yCoord_vec  = yCoordGeo_vec;
    xCoord_high = xCoordGeo_vec;
    yCoord_high = yCoordGeo_vec;
end

if order_G > 1
    [quadG_high_conn, xCoordG_high, yCoordG_high] = quad4_to_quadHigh(connectivityMatrixGeo_mat, xCoordGeo_vec, yCoordGeo_vec, order_G);
    connectivityMatrixGeo_mat = quadG_high_conn;
    xCoordGeo_vec  = xCoordG_high;
    yCoordGeo_vec  = yCoordG_high;

end

%% Precompute geometric data for all elements 
numElements = numElementsX * numElementsY;
elementData = cell(numElements, 1);

for eleNum = 1:numElements

    elementData{eleNum}.GNodes           = connectivityMatrixGeo_mat(eleNum,:);
    elementData{eleNum}.xNodesGVals_vec  = xCoordGeo_vec(elementData{eleNum}.GNodes);
    elementData{eleNum}.yNodesGVals_vec  = yCoordGeo_vec(elementData{eleNum}.GNodes);

    elementData{eleNum}.Nodes            = connectivityMatrix_mat(eleNum,:);
    elementData{eleNum}.xNodesVals_vec   = xCoord_vec(elementData{eleNum}.Nodes);
    elementData{eleNum}.yNodesVals_vec   = yCoord_vec(elementData{eleNum}.Nodes);
    
    [weights_pages,N_row_points_pages,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,ElementGeoCoords_mat] = getGaussRelated(elementData{eleNum}.xNodesGVals_vec,elementData{eleNum}.yNodesGVals_vec,K,numGaussPoints,K_G,numGeometryGaussPoints,order,order_G);

    elementData{eleNum}.weights    = weights_pages;
    elementData{eleNum}.N_row      = N_row_points_pages;
    elementData{eleNum}.J_det      = J_det_points_elev;
    elementData{eleNum}.dNdX       = N_diff_PhysCoords_points_rows_pages(1, :, :);
    elementData{eleNum}.dNdY       = N_diff_PhysCoords_points_rows_pages(2, :, :);
    elementData{eleNum}.physCoords = ElementGeoCoords_mat;
end

%% Compute Global Matrices
numElements = size(connectivityMatrix_mat,1);
totNumNodes = length(xCoord_vec);
pK_mat      = zeros(totNumNodes,totNumNodes);
uMLV_col    = zeros(totNumNodes,1);
vMLV_col    = zeros(totNumNodes,1);

for eleNum = 1:size(connectivityMatrix_mat,1)

    elementNodes_vec = elementData{eleNum}.Nodes;
    xNodesVals_vec   = elementData{eleNum}.xNodesVals_vec;
    yNodesVals_vec   = elementData{eleNum}.yNodesVals_vec;

    elementNodesGeo_vec = elementData{eleNum}.GNodes;
    xNodesValsGeo_vec   = elementData{eleNum}.xNodesGVals_vec;
    yNodesValsGeo_vec   = elementData{eleNum}.yNodesGVals_vec;

    [pk_mat,umlv_col,vmlv_col] = computeLocalNSP(elementData{eleNum});

    pK_mat(elementNodes_vec,elementNodes_vec) = pK_mat(elementNodes_vec,elementNodes_vec) + pk_mat;
    uMLV_col(elementNodes_vec)                = uMLV_col(elementNodes_vec)                + umlv_col';
    vMLV_col(elementNodes_vec)                = vMLV_col(elementNodes_vec)                + vmlv_col';
end

%% Boundary Conditions locations
pBCCondition_vec       = (xCoord_vec == lengthX/2) & (yCoord_vec == 0);
lowerWall_vec          = (yCoord_vec == 0);
upperWallNoCorners_vec = (yCoord_vec == lengthY) & (xCoord_vec ~= lengthX) & (xCoord_vec ~= 0);
rightWallNoCorner_vec  = (xCoord_vec == lengthX) & (yCoord_vec ~= lengthY);
leftWallNoCorner_vec   = (xCoord_vec == 0)       & (yCoord_vec ~= lengthY);
upperLeftCorner_vec    = (xCoord_vec == 0)       & (yCoord_vec == lengthY);
upperRightCorner_vec   = (xCoord_vec == lengthX) & (yCoord_vec == lengthY);
xEqualHalf             =  xCoord_vec == lengthX/2;

%% Prepare Pressure Global Matrix
ind_vec                           = 1:length(xCoord_vec);
pBCNodes_vec                      = ind_vec(pBCCondition_vec);
pK_mat(pBCNodes_vec,:)            = 0;
pK_mat(pBCNodes_vec,pBCNodes_vec) = eye(length(pBCNodes_vec));

%% Initial Conditions
un_col = zeros(totNumNodes,1);
vn_col = zeros(totNumNodes,1);

%% time Loop
itr  = 1;
time = 0;
pK_mat  = decomposition(pK_mat);   % this line to enhance performance

while(itr < totitrs)
    pressureSolution_col = computePressure(elementData,connectivityMatrix_mat,pBCNodes_vec,pK_mat,totNumNodes,un_col,vn_col,epslon);

    uRHS_col    = zeros(totNumNodes,1);
    vRHS_col    = zeros(totNumNodes,1);
    ReLocal_col = zeros(numElements,1);
    CFL_col     = zeros(numElements,1);
    xCenter_col = zeros(numElements,1);
    yCenter_col = zeros(numElements,1);

    for eleNum = 1:numElements

        elementNodes_vec = elementData{eleNum}.Nodes;

        unNodes_col                = un_col(elementNodes_vec,1);
        vnNodes_col                = vn_col(elementNodes_vec,1);
        pressureSolutionNodes_col  = pressureSolution_col(elementNodes_vec,1);

        [urhs_col,vrhs_col,ReLocal,cfl,xCenter,yCenter] = computeLocalNS(elementData{eleNum},unNodes_col,vnNodes_col,pressureSolutionNodes_col,Re,timeStep);

        uRHS_col(elementNodes_vec)  = uRHS_col(elementNodes_vec)   + urhs_col;
        vRHS_col(elementNodes_vec)  = vRHS_col(elementNodes_vec)   + vrhs_col;

        ReLocal_col(eleNum)  = ReLocal;
        CFL_col(eleNum)      = cfl;
        xCenter_col(eleNum)  = xCenter;
        yCenter_col(eleNum)  = yCenter;
    end

    %% u-velocity Component Boundary Condition
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

    %% v-velocity Component Boundary Condition
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

    %% u Solution
    uSolution_col = uRHS_col./uMLV_col;

    %% v Solution
    vSolution_col = vRHS_col./vMLV_col;

    %% Reassigning initial Conditions
    un_col  = uSolution_col;
    vn_col  = vSolution_col;

    if mod(itr, vtkOutputInterval) == 0
        writeVTKatInterval(outputFolder,itr,xCoord_vec,yCoord_vec,un_col,vn_col,connectivityMatrix_mat,pressureSolution_col,time)
    end

    %% Stepping
    time = time + timeStep;
    itr  = itr + 1;
    disp(itr)
end
toc
writeFinalVTK(xCoord_vec,yCoord_vec,un_col,vn_col,connectivityMatrix_mat,pressureSolution_col,time,outputFolder)

%% Save results
resultsFileName = sprintf('cavity_resultsRe%d.mat', Re);
save(resultsFileName, 'uSolution_col', 'vSolution_col', 'pressureSolution_col', ...
    'xCoord_vec', 'yCoord_vec', 'connectivityMatrix_mat', ...
    'xCoordGeo_vec', 'yCoordGeo_vec', 'connectivityMatrixGeo_mat');
