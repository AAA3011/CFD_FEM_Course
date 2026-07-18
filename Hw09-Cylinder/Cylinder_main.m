% main - Main simulation driver.
%
%
% Executes transient FEM solver with time-stepping, reads geometry, precomputes
% element data, assembles global matrices, solves the pressure and velocity
% equations, computes forces, and writes results.
%
% Inputs: global parameters set in script (epslon, Re, timeStep, vtkOutputInterval,
% outputFolder).
% Outputs: `results.mat` and VTK files (`*.vtu`) written to the output folder.
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
epslon    = 0.05;
Re        = 100;
timeStep  = 0.005;

vtkOutputInterval = 10;
outputFolder      = 'Cylinder100Re';
makeVTKdir(outputFolder)

%% Reading Geometry and BC
[connectivityMatrix_mat,x_col,y_col,top_vec,outlet_vec,inlet_vec,cylinderWall_vec,bottom_vec] = readGeoAndBC();

%% Precompute geometric data for all elements
numGaussPoints = 5;
n_ElementType  = 1;

numElements    = size(connectivityMatrix_mat,1);
totNumNodes    = length(x_col);
elementData    = cell(numElements, 1);

for eleNum = 1:numElements

    elementData{eleNum}.Nodes           = connectivityMatrix_mat(eleNum,:);
    elementData{eleNum}.xNodesVals_vec  = x_col(elementData{eleNum}.Nodes);
    elementData{eleNum}.yNodesVals_vec  = y_col(elementData{eleNum}.Nodes);

    [weights_pages,N_row_points_pages,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,ElementPhysCoords_mat] = getGaussRelated(elementData{eleNum}.xNodesVals_vec,elementData{eleNum}.yNodesVals_vec,n_ElementType,numGaussPoints);

    elementData{eleNum}.weights    = weights_pages;
    elementData{eleNum}.N_row      = N_row_points_pages;
    elementData{eleNum}.J_det      = J_det_points_elev;
    elementData{eleNum}.dNdX       = N_diff_PhysCoords_points_rows_pages(1, :, :);
    elementData{eleNum}.dNdY       = N_diff_PhysCoords_points_rows_pages(2, :, :);
    elementData{eleNum}.physCoords = ElementPhysCoords_mat;

end

%% Compute Global Matrices
% pK_mat      = zeros(totNumNodes,totNumNodes);
% uMLV_col    = zeros(totNumNodes,1);
% vMLV_col    = zeros(totNumNodes,1);
% MLV_col     = zeros(totNumNodes,1);
% 
% for eleNum = 1:size(connectivityMatrix_mat,1)
% 
%     elementNodes_vec = elementData{eleNum}.Nodes;
% 
%     [pk_mat,umlv_col,vmlv_col,mlv_col] = computeLocalNSP(elementData{eleNum});
% 
%     pK_mat(elementNodes_vec,elementNodes_vec) = pK_mat(elementNodes_vec,elementNodes_vec) + pk_mat;
%     uMLV_col(elementNodes_vec)                = uMLV_col(elementNodes_vec)                + umlv_col';
%     vMLV_col(elementNodes_vec)                = vMLV_col(elementNodes_vec)                + vmlv_col';
%     MLV_col(elementNodes_vec)                 = MLV_col(elementNodes_vec)                 + mlv_col';
% end

%% Compute Global Matrices
numNodesPerEle = size(connectivityMatrix_mat, 2); 
totalEntries   = numElements * (numNodesPerEle^2);

% Preallocate triplet vectors for the sparse matrix
I = zeros(totalEntries, 1);
J = zeros(totalEntries, 1);
V = zeros(totalEntries, 1);

uMLV_col    = zeros(totNumNodes,1);
vMLV_col    = zeros(totNumNodes,1);
MLV_col     = zeros(totNumNodes,1);

cnt = 0;
for eleNum = 1:numElements
    elementNodes_vec = elementData{eleNum}.Nodes;
    [pk_mat,umlv_col,vmlv_col,mlv_col] = computeLocalNSP(elementData{eleNum});

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
    MLV_col(elementNodes_vec)   = MLV_col(elementNodes_vec)  + mlv_col';
end

% Create the Global Sparse Matrix directly without high memory overhead
pK_mat = sparse(I, J, V, totNumNodes, totNumNodes);

%% Prepare Pressure Global Matrix
[~, localIdx]  = min(abs(y_col(outlet_vec)));
pBCNodes_col = outlet_vec(localIdx);

pK_mat(pBCNodes_col,:)              = 0;
pK_mat(pBCNodes_col,pBCNodes_col) = eye(length(pBCNodes_col));


%% Initial Conditions and Time
time   = 0;
un_col = ones(totNumNodes,1);
vn_col = zeros(totNumNodes,1);

%% Error parameters
itr     = 1;
totitrs = 50000;
Fx_col = zeros(totitrs,1);
Fy_col = zeros(totitrs,1);

%% Decomposition
pK_mat  = decomposition(pK_mat);   % this line to enhance performance

while(itr < totitrs)
    pressureSolution_col = computePressure(elementData,numElements,pBCNodes_col,pK_mat,totNumNodes,un_col,vn_col,epslon);

    %% Initialize Matrices
    uRHS_col    = zeros(totNumNodes,1);
    vRHS_col    = zeros(totNumNodes,1);
    uxRHS_col   = zeros(totNumNodes,1);
    vxRHS_col   = zeros(totNumNodes,1);
    uyRHS_col   = zeros(totNumNodes,1);
    vyRHS_col   = zeros(totNumNodes,1);
    ReLocal_col = zeros(numElements,1);
    CFL_col     = zeros(numElements,1);
    xCenter_col = zeros(numElements,1);
    yCenter_col = zeros(numElements,1);

    for elementNumber = 1:numElements

        elementNodes_vec           = elementData{elementNumber}.Nodes;

        unNodes_col                = un_col(elementNodes_vec,1);
        vnNodes_col                = vn_col(elementNodes_vec,1);
        pressureSolutionNodes_col  = pressureSolution_col(elementNodes_vec,1);

        [urhs_col,vrhs_col,ReLocal,cfl,xCenter,yCenter,uxrhs_col,uyrhs_col,vxrhs_col,vyrhs_col] = computeLocalNS(elementData{elementNumber},unNodes_col,vnNodes_col,pressureSolutionNodes_col,Re,timeStep);

        uRHS_col(elementNodes_vec)  = uRHS_col(elementNodes_vec)  + urhs_col;
        vRHS_col(elementNodes_vec)  = vRHS_col(elementNodes_vec)  + vrhs_col;
        uxRHS_col(elementNodes_vec) = uxRHS_col(elementNodes_vec) + uxrhs_col';
        vxRHS_col(elementNodes_vec) = vxRHS_col(elementNodes_vec) + vxrhs_col';
        uyRHS_col(elementNodes_vec) = uyRHS_col(elementNodes_vec) + uyrhs_col';
        vyRHS_col(elementNodes_vec) = vyRHS_col(elementNodes_vec) + vyrhs_col';
        ReLocal_col(elementNumber)  = ReLocal;
        CFL_col(elementNumber)      = cfl;
        xCenter_col(elementNumber)  = xCenter;
        yCenter_col(elementNumber)  = yCenter;

    end

    %% u-velocity Component Boundary Condition
    uRHS_col(inlet_vec)        = 1;
    uRHS_col(cylinderWall_vec) = 0;

    uMLV_col(inlet_vec)        = 1;
    uMLV_col(cylinderWall_vec) = 1;

    %% v-velocity Component Boundary Condition
    vRHS_col(inlet_vec)        = 0;
    vRHS_col(top_vec)          = 0;
    vRHS_col(bottom_vec)       = 0;
    vRHS_col(cylinderWall_vec) = 0;

    vMLV_col(inlet_vec)        = 1;
    vMLV_col(top_vec)          = 1;
    vMLV_col(bottom_vec)       = 1;
    vMLV_col(cylinderWall_vec) = 1;

    %% Solve for u velocity
    uSolution_col = uRHS_col./uMLV_col;

    %% Solve for v velocity
    vSolution_col = vRHS_col./vMLV_col;
    %% Solve for ux
    uxSolution_col= uxRHS_col./MLV_col;

    %% Solve for uy
    uySolution_col= uyRHS_col./MLV_col;

    %% Solve for vx
    vxSolution_col= vxRHS_col./MLV_col;

    %% Solve for vy
    vySolution_col= vyRHS_col./MLV_col;

    %% Forces
    [Fx,Fy]     = compForces(uxSolution_col,uySolution_col,vxSolution_col,vySolution_col,pressureSolution_col,x_col,y_col,cylinderWall_vec,Re);
    Fx_col(itr) = Fx;
    Fy_col(itr) = Fy;

    %% Write VTK file at specified intervals
    if mod(itr, vtkOutputInterval) == 0
        writeVTKatInterval(outputFolder,itr,x_col,y_col,un_col,vn_col,connectivityMatrix_mat,pressureSolution_col,time)
    end

    %% Reassigning initial Conditions
    un_col  = uSolution_col;
    vn_col  = vSolution_col;

    %% Stepping
    time = time + timeStep;
    %% iteration
    itr = itr + 1;
    disp(itr)
end
toc

% Save results
outputFolder = 'Cylinder100Re';
save(fullfile(outputFolder, 'results.mat'), 'Fx_col', 'Fy_col', 'un_col', 'vn_col', 'pressureSolution_col', 'time');

%% Final VTK file
writeFinalVTK(x_col,y_col,un_col,vn_col,connectivityMatrix_mat,pressureSolution_col,time,outputFolder)


