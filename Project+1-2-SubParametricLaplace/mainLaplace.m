% mainLaplace - Main simulation driver.
%
% FILE: mainLaplace.m
% DESCRIPTION:
% Driver script that sets up a structured mesh, constructs high-order
% quadrilateral elements, precomputes element data, assembles the global
% stiffness matrix for a Laplace problem, applies boundary conditions,
% solves the linear system, computes error norms, and plots results.
%
% INPUTS (script parameters inside file):
% numGaussPoints, order, order_G, numElementsX, numElementsY, lengthX/Y
%
% OUTPUTS:
% Figures saved to ./figures and computed variables in workspace
%
% Syntax: mainLaplace()
%
% Inputs: none
%
% Outputs: none
%

clc;clearvars;close all;

%% Element informations
numGaussPoints         = 10;
numGeometryGaussPoints = numGaussPoints;

order   = 5;
K       = (order+1)^2;

order_G = 1;
K_G     = (order_G+1)^2;

%% Connectivity Matrix and Grid Generation
numElementsX = 14;
numElementsY = 14;
numNodesX    = numElementsX + 1;
numNodesY    = numElementsY + 1;

lengthX      = 1;
lengthY      = 1;
xMin         = 0;
yMin         = 0;

[connectivityMatrixGeo_mat,xCoordGeo_vec,yCoordGeo_vec] = structuredMesh(numElementsX,numElementsY,lengthX,lengthY,xMin,yMin);

%% Converting to higher order mesh
[quad_high_conn, xCoord_high, yCoord_high] = quad4_to_quadHigh(connectivityMatrixGeo_mat, xCoordGeo_vec, yCoordGeo_vec, order);
connectivityMatrix_mat = quad_high_conn;
xCoord_vec = xCoord_high;
yCoord_vec = yCoord_high;

if order_G > 1
    [quad_high_conn, xCoord_high, yCoord_high] = quad4_to_quadHigh(connectivityMatrixGeo_mat, xCoordGeo_vec, yCoordGeo_vec, order);
    connectivityMatrixGeo_mat = quad_high_conn;
    xCoordGeo_vec = xCoord_high;
    yCoordGeo_vec = yCoord_high;
end

%% Precompute geometric data for all elements 
numElements = numElementsX * numElementsY;
elementData = cell(numElements, 1);

for eleNum = 1:numElements
    elementData{eleNum}.GNodes             = connectivityMatrixGeo_mat(eleNum,:);
    elementData{eleNum}.xNodesValsGeo_vec  = xCoordGeo_vec(elementData{eleNum}.GNodes);
    elementData{eleNum}.yNodesValsGeo_vec  = yCoordGeo_vec(elementData{eleNum}.GNodes);

    elementData{eleNum}.Nodes           = connectivityMatrix_mat(eleNum,:);
    elementData{eleNum}.xNodesVals_vec  = xCoord_vec(elementData{eleNum}.Nodes);
    elementData{eleNum}.yNodesVals_vec  = yCoord_vec(elementData{eleNum}.Nodes);

    [weights_pages,N_row_points_pages,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,ElementPhysCoords_mat] = getGaussRelated(elementData{eleNum}.xNodesValsGeo_vec,elementData{eleNum}.yNodesValsGeo_vec,K,numGaussPoints,K_G,numGeometryGaussPoints,order,order_G);


    elementData{eleNum}.weights    = weights_pages;
    elementData{eleNum}.N_row      = N_row_points_pages;
    elementData{eleNum}.J_det      = J_det_points_elev;
    elementData{eleNum}.dNdX       = N_diff_PhysCoords_points_rows_pages(1, :, :);
    elementData{eleNum}.dNdY       = N_diff_PhysCoords_points_rows_pages(2, :, :);
    elementData{eleNum}.physCoords = ElementPhysCoords_mat;
end

%% local, Global Matrices
totNumNodes      = length(xCoord_vec);
totNumElements   = size(connectivityMatrix_mat,1);
globalMatrix_mat = zeros(totNumNodes);

for eleNum = 1:totNumElements

    elementNodes_vec = elementData{eleNum}.Nodes;

    elementNodesGeo_vec = elementData{eleNum}.GNodes;
    xNodesValsGeo_vec   = elementData{eleNum}.xNodesValsGeo_vec;
    yNodesValsGeo_vec   = elementData{eleNum}.yNodesValsGeo_vec;

    localMatrix_mat = calcLocal(elementData{eleNum});

    globalMatrix_mat(elementNodes_vec,elementNodes_vec) = globalMatrix_mat(elementNodes_vec,elementNodes_vec) + localMatrix_mat;

end

%% Boundary Conditions Locations
ind_vec    = 1:totNumNodes;
upperWall  = ind_vec(yCoord_vec == 1);
bottomWall = ind_vec(yCoord_vec == 0);
rightWall  = ind_vec(xCoord_vec == 1);
leftWall   = ind_vec(xCoord_vec == 0);

%% Forcing Boundary Conditions
uRHS_cvec = zeros(totNumNodes,1);
uRHS_cvec(bottomWall) = sin(pi * xCoord_vec(bottomWall));

globalMatrix_mat(bottomWall,:)          = 0;
globalMatrix_mat(bottomWall,bottomWall) = eye(length(bottomWall));
globalMatrix_mat(upperWall,:)           = 0;
globalMatrix_mat(upperWall,upperWall)   = eye(length(upperWall));
globalMatrix_mat(rightWall,:)           = 0;
globalMatrix_mat(rightWall,rightWall)   = eye(length(rightWall));
globalMatrix_mat(leftWall,:)            = 0;
globalMatrix_mat(leftWall,leftWall)     = eye(length(leftWall));

%% FDM Solution
uSolution_cvec = globalMatrix_mat \ uRHS_cvec;

%% Exact Solution
uExact_cvec  = (cosh(pi*yCoord_vec) - coth(pi)*sinh(pi*yCoord_vec)).* sin(pi*xCoord_vec);
uxExact_cvec = pi .* ( cosh(pi*yCoord_vec) - coth(pi).*sinh(pi*yCoord_vec) ) .* cos(pi*xCoord_vec);
uyExact_cvec = pi*( sinh(pi*yCoord_vec) - coth(pi)*cosh(pi*yCoord_vec) ).* sin(pi*xCoord_vec);

%% Errors Calculations
[L2_norm,H1_norm] =  calcNorms(uSolution_cvec,totNumElements,elementData);

%% Plotting Code
fsLabel = 14;
fsTitle = 16;
fsLegend = 12;
fsAxes = 12;

outDir = fullfile(pwd, 'figures');
if ~exist(outDir, 'dir')
    [ok, msg] = mkdir(outDir);
    if ~ok
        error('Failed to create folder %s: %s', outDir, msg);
    end
end

% Solution Plotting: FEM vs Exact, side by side
tri = delaunay(xCoord_vec, yCoord_vec);
figFEMExact = figure('Position',[100 100 1200 500]);

subplot(1,2,1)
trisurf(tri, xCoord_vec, yCoord_vec, uSolution_cvec)
shading interp
view(2)
colorbar
axis equal tight
title(['FEM solution using ', num2str(totNumElements),' elements Q', num2str(K_G),' geometry Q', num2str(K), ' solution'], ...
    'Interpreter','latex','FontSize',fsTitle)
set(gca,'FontSize',fsAxes)

subplot(1,2,2)
trisurf(tri, xCoord_vec, yCoord_vec, uExact_cvec)
view(2); shading interp; colorbar
axis equal tight
title('Exact solution','Interpreter','latex','FontSize',fsTitle)
set(gca,'FontSize',fsAxes)

print(figFEMExact, fullfile(outDir,'FEM_vs_Exact_solution'), '-dsvg')

% Solution plotting 3D
fig3D = figure;
trisurf(delaunay(xCoord_vec,yCoord_vec),xCoord_vec, yCoord_vec, uSolution_cvec)
shading interp
colorbar
axis tight
xlabel('x','FontSize',fsLabel); ylabel('y','FontSize',fsLabel); zlabel('u','FontSize',fsLabel)
title('FEM solution','Interpreter','latex','FontSize',fsTitle)
set(gca,'FontSize',fsAxes)
print(fig3D, fullfile(outDir,'FEM_solution_3D'), '-dsvg')

% Mesh Plotting
figMesh = figure('Position',[100 100 1400 700]*0.85);
subplot(1,2,1)
plotQuadHigherMesh(connectivityMatrixGeo_mat, xCoordGeo_vec, yCoordGeo_vec)
axis equal;
xlabel('x/d','Interpreter','latex','FontSize',fsLabel);
ylabel('y/d','Interpreter','latex','FontSize',fsLabel);
title(['Geometry mesh Q', num2str(K_G)],'Interpreter','latex','FontSize',fsTitle);
grid on;
set(gca,'FontSize',fsAxes)

subplot(1,2,2)
plotQuadHigherMesh(quad_high_conn, xCoord_high, yCoord_high)
axis equal;
xlabel('x/d','Interpreter','latex','FontSize',fsLabel);
ylabel('y/d','Interpreter','latex','FontSize',fsLabel);
title(['Solution mesh Q', num2str(K)],'Interpreter','latex','FontSize',fsTitle);
grid on;
set(gca,'FontSize',fsAxes)
% print(figMesh, fullfile(outDir,'mesh_plot'), '-dsvg')
print(figMesh, fullfile(outDir,'mesh_plot'), '-dpng', '-r300')

% Solution along x = 0.5
figLine = figure;
xTarget = 0.5;
tolLine = 1e-10;
lineIdx = find(abs(xCoord_vec - xTarget) < tolLine);
if isempty(lineIdx)
    warning('No nodes found near x = 0.5. Increase tolLine or check mesh.');
else
    [yLine, sortIdx] = sort(yCoord_vec(lineIdx));
    uSolLine = uSolution_cvec(lineIdx(sortIdx));
    uExLine = uExact_cvec(lineIdx(sortIdx));
    plot(yLine, uSolLine, 'LineWidth', 2, 'DisplayName', '$$u_{sol}$$');
    hold on;
    plot(yLine, uExLine,'o', 'LineWidth', 1.2, 'DisplayName', '$$u_{exact}$$');
    xlabel('y','Interpreter','latex','FontSize',fsLabel);
    ylabel('u(x=0.5,y)','Interpreter','latex','FontSize',fsLabel);
    title('u values along x = 0.5','Interpreter','latex','FontSize',fsTitle);
    legend('Interpreter','latex','FontSize',fsLegend);
    grid on;
    set(gca,'FontSize',fsAxes)
    print(figLine, fullfile(outDir,'solution_along_x_0p5'), '-dsvg')
end
