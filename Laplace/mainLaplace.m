clc;clearvars;close all;

%% Element informations
numGaussPoints         = 10;
numGeometryGaussPoints = numGaussPoints;

order   = 4;
K       = (order+1)^2;

order_G = 1;
K_G     = (order_G+1)^2;

%% Connectivity Matrix and Grid Generation
numElementsX = 2;
numElementsY = 2;
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

%% local, Global Matrices
totNumNodes      = length(xCoord_vec);
totNumElements   = size(connectivityMatrix_mat,1);
globalMatrix_mat = zeros(totNumNodes);

for eleNum = 1:totNumElements

    elementNodes_vec = connectivityMatrix_mat(eleNum,:);

    elementNodesGeo_vec = connectivityMatrixGeo_mat(eleNum,:);
    xNodesValsGeo_vec   = xCoordGeo_vec(elementNodesGeo_vec);
    yNodesValsGeo_vec   = yCoordGeo_vec(elementNodesGeo_vec);

    localMatrix_mat = calcLocal(numGaussPoints,K,K_G,numGeometryGaussPoints,xNodesValsGeo_vec,yNodesValsGeo_vec,order,order_G);

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

globalMatrix_mat(bottomWall,:) = 0;
globalMatrix_mat(bottomWall,bottomWall) = eye(length(bottomWall));
globalMatrix_mat(upperWall,:) = 0;
globalMatrix_mat(upperWall,upperWall) = eye(length(upperWall));
globalMatrix_mat(rightWall,:) = 0;
globalMatrix_mat(rightWall,rightWall) = eye(length(rightWall));
globalMatrix_mat(leftWall,:) = 0;
globalMatrix_mat(leftWall,leftWall) = eye(length(leftWall));

%% FDM Solution
uSolution_cvec = globalMatrix_mat \ uRHS_cvec;

%% Exact Solution
uExact_cvec  = (cosh(pi*yCoord_vec) - coth(pi)*sinh(pi*yCoord_vec)).* sin(pi*xCoord_vec);
uxExact_cvec = pi .* ( cosh(pi*yCoord_vec) - coth(pi).*sinh(pi*yCoord_vec) ) .* cos(pi*xCoord_vec);
uyExact_cvec = pi*( sinh(pi*yCoord_vec) - coth(pi)*cosh(pi*yCoord_vec) ).* sin(pi*xCoord_vec);

%% Errors Calculations
[L2_norm,H1_norm] =  calcNorms(uSolution_cvec,connectivityMatrix_mat,connectivityMatrixGeo_mat,totNumElements,xCoordGeo_vec,yCoordGeo_vec,K,K_G,order,order_G,numGeometryGaussPoints,numGaussPoints,xCoord_vec,yCoord_vec);

%% Plotting Code
fsLabel = 14;
fsTitle = 16;
fsLegend = 12;
fsAxes = 12;
% Solution Plotting
figure;
tri = delaunay(xCoord_vec, yCoord_vec);
trisurf(tri, xCoord_vec, yCoord_vec, uSolution_cvec)
shading interp
view(2)
colorbar
axis equal tight
title(['FEM solution using ', num2str(totNumElements),' elements Q', num2str(K_G),' geometry Q', num2str(K), ' solution'], ...
    'Interpreter','latex','FontSize',fsTitle)
set(gca,'FontSize',fsAxes)

figure;
trisurf(tri, xCoord_vec, yCoord_vec, uExact_cvec)
view(2); shading interp; colorbar
axis equal tight
title('Exact solution','Interpreter','latex','FontSize',fsTitle)
set(gca,'FontSize',fsAxes)

% Solution plotting 3D
figure
trisurf(delaunay(xCoord_vec,yCoord_vec),xCoord_vec, yCoord_vec, uSolution_cvec)
shading interp
colorbar
axis tight
xlabel('x','FontSize',fsLabel); ylabel('y','FontSize',fsLabel); zlabel('u','FontSize',fsLabel)
title('FEM solution','Interpreter','latex','FontSize',fsTitle)
set(gca,'FontSize',fsAxes)

% Mesh Plotting
figure;
plotQuadHigherMesh(quad_high_conn,xCoord_high, yCoord_high)
axis equal;
xlabel('x/d','Interpreter','latex','FontSize',fsLabel);
ylabel('y/d','Interpreter','latex','FontSize',fsLabel);
title(['Mesh for Q' , num2str(K), ' Element'] ,'Interpreter','latex','FontSize',fsTitle);
grid on;
set(gca,'FontSize',fsAxes)

% Solution along x = 0.5
figure;
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
end

