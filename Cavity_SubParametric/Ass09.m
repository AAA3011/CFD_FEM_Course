clc;clearvars;close all;
tic
%% Element informations
numGaussPoints = 10;
numGeometryGaussPoints = 10;
K     = 16;
order = 3;
K_G   = 4;

%% Global Parameters
Re        = 100;
epslon    = 0.01;
timeStep  = 0.005;
totitrs   = 1000;

vtkOutputInterval = 10;
outputFolder      = 'CavityResutts100ReHigher';
makeVTKdir(outputFolder)

%% Connectivity Matrix and Grid Generation
numElementsX = 1;
numElementsY = 1;
% numNodesX    = numElementsX + 1;
% numNodesY    = numElementsY + 1;
% totNumNodes  = numNodesX * numNodesY;
lengthX      = 1;
lengthY      = 1;
xMin         = 0;
yMin         = 0;

[connectivityMatrixGeo_mat,xCoordGeo_vec,yCoordGeo_vec] = structuredMesh(numElementsX,numElementsY,lengthX,lengthY,xMin,yMin);
% [quad9_conn, xCoord9, yCoord9]                          = quad4_to_quad9(connectivityMatrixGeo_mat, xCoordGeo_vec, yCoordGeo_vec);
% connectivityMatrix_mat = quad9_conn;
% xCoord_vec = xCoord9';
% yCoord_vec = yCoord9';

[quad_high_conn, xCoord_high, yCoord_high] = quad4_to_quadHigh(connectivityMatrixGeo_mat, xCoordGeo_vec, yCoordGeo_vec, order);

% [quad_high_conn, xCoord_high, yCoord_high]   = quad4_to_quad9(connectivityMatrixGeo_mat, xCoordGeo_vec, yCoordGeo_vec);
connectivityMatrix_mat = quad_high_conn;
xCoord_vec = xCoord_high';
yCoord_vec = yCoord_high';

% connectivityMatrixGeo_mat = quad_high_conn;
% xCoordGeo_vec = xCoord_high';
% yCoordGeo_vec = yCoord_high';
plot_quad9_mesh(quad_high_conn,xCoord_high, yCoord_high)

totNumNodes = length(xCoord_vec);
%% Compute Global Matrices
numElements  = size(connectivityMatrix_mat,1);
pK_mat       = zeros(totNumNodes,totNumNodes);
uMLV_cvec    = zeros(totNumNodes,1);
vMLV_cvec    = zeros(totNumNodes,1);
for eleNum = 1:size(connectivityMatrix_mat,1)

    elementNodes_vec = connectivityMatrix_mat(eleNum,:);
    xNodesVals_vec   = xCoord_vec(elementNodes_vec);
    yNodesVals_vec   = yCoord_vec(elementNodes_vec);

    elementNodesGeo_vec = connectivityMatrixGeo_mat(eleNum,:);
    xNodesValsGeo_vec   = xCoord_vec(elementNodesGeo_vec);
    yNodesValsGeo_vec   = yCoord_vec(elementNodesGeo_vec);

    [pk_mat,umlv_cvec,vmlv_cvec] = computeLocalNSP(numGaussPoints,K,xNodesVals_vec,yNodesVals_vec,K_G,numGeometryGaussPoints,xNodesValsGeo_vec,yNodesValsGeo_vec);

    pK_mat(elementNodes_vec,elementNodes_vec) = pK_mat(elementNodes_vec,elementNodes_vec) + pk_mat;
    uMLV_cvec(elementNodes_vec)  = uMLV_cvec(elementNodes_vec)   + umlv_cvec';
    vMLV_cvec(elementNodes_vec)  = vMLV_cvec(elementNodes_vec)   + vmlv_cvec';
end

%% Boundary Conditions locations
lengthX = 1;
lengthY = 1;
ind_vec = 1:length(xCoord_vec);

pBCCondition_vec       = (xCoord_vec == lengthX/2) & (yCoord_vec == 0);
lowerWall_vec          = (yCoord_vec == 0);
upperWallNoCorners_vec = (yCoord_vec == lengthY) & (xCoord_vec ~= lengthX) & (xCoord_vec ~= 0);
rightWallNoCorner_vec  = (xCoord_vec == lengthX) & (yCoord_vec ~= lengthY);
leftWallNoCorner_vec   = (xCoord_vec == 0)       & (yCoord_vec ~= lengthY);
upperLeftCorner_vec    = (xCoord_vec == 0)       & (yCoord_vec == lengthY);
upperRightCorner_vec   = (xCoord_vec == lengthX) & (yCoord_vec == lengthY);
xEqualHalf             =  xCoord_vec == 0.5;

%% Prepare Pressure Global Matrix
pBCNodes_vec                      = ind_vec(pBCCondition_vec);
pK_mat(pBCNodes_vec,:)            = 0;
pK_mat(pBCNodes_vec,pBCNodes_vec) = eye(length(pBCNodes_vec));

%% Initial Conditions
time      = 0;
un_cvec = zeros(totNumNodes,1);
vn_cvec = zeros(totNumNodes,1);

%% Error parameters
itr = 1;
while(itr < totitrs)

    pressureSolution_cvec = computePressure(numGaussPoints,K,connectivityMatrix_mat,xCoord_vec,yCoord_vec,pBCNodes_vec,pK_mat,totNumNodes,un_cvec,vn_cvec,epslon,K_G,numGeometryGaussPoints,connectivityMatrixGeo_mat);

    uRHS_cvec    = zeros(totNumNodes,1);
    vRHS_cvec    = zeros(totNumNodes,1);
    ReLocal_cvec = zeros(numElements,1);
    CFL_cvec     = zeros(numElements,1);
    xCenter_cvec = zeros(numElements,1);
    yCenter_cvec = zeros(numElements,1);

    for elementNumber = 1:numElements

        elementNodesGeo_vec = connectivityMatrixGeo_mat(elementNumber,:);
        xNodesValsGeo_vec   = xCoord_vec(elementNodesGeo_vec);
        yNodesValsGeo_vec   = yCoord_vec(elementNodesGeo_vec);

        elementNodes_vec = connectivityMatrix_mat(elementNumber,:);
        xNodesVals_vec   = xCoord_vec(elementNodes_vec);
        yNodesVals_vec   = yCoord_vec(elementNodes_vec);
        unNodes_cvec     = un_cvec(elementNodes_vec,1);
        vnNodes_cvec     = vn_cvec(elementNodes_vec,1);
        pressureSolutionNodes_cvec  = pressureSolution_cvec(elementNodes_vec,1);
        [urhs_cvec,vrhs_cvec,ReLocal,cfl,xCenter,yCenter] = computeLocalNS(numGaussPoints,K,xNodesVals_vec,yNodesVals_vec,unNodes_cvec,vnNodes_cvec,pressureSolutionNodes_cvec,Re,timeStep,K_G,numGeometryGaussPoints,xNodesValsGeo_vec,yNodesValsGeo_vec);


        uRHS_cvec(elementNodes_vec)  = uRHS_cvec(elementNodes_vec)   + urhs_cvec;
        vRHS_cvec(elementNodes_vec)  = vRHS_cvec(elementNodes_vec)   + vrhs_cvec;
        ReLocal_cvec(elementNumber)  = ReLocal;
        CFL_cvec(elementNumber)      = cfl;
        xCenter_cvec(elementNumber)  = xCenter;
        yCenter_cvec(elementNumber)  = yCenter;

    end

    %% u-velocity Component Boundary Condition
    uRHS_cvec(ind_vec(lowerWall_vec))          = 0;
    uRHS_cvec(ind_vec(upperWallNoCorners_vec)) = 1;
    uRHS_cvec(ind_vec(rightWallNoCorner_vec))  = 0;
    uRHS_cvec(ind_vec(leftWallNoCorner_vec))   = 0;
    uRHS_cvec(ind_vec(upperLeftCorner_vec))    = 1;
    uRHS_cvec(ind_vec(upperRightCorner_vec))   = 1;

    uMLV_cvec(ind_vec(lowerWall_vec))          = 1;
    uMLV_cvec(ind_vec(upperWallNoCorners_vec)) = 1;
    uMLV_cvec(ind_vec(rightWallNoCorner_vec))  = 1;
    uMLV_cvec(ind_vec(leftWallNoCorner_vec))   = 1;
    uMLV_cvec(ind_vec(upperLeftCorner_vec))    = 1;
    uMLV_cvec(ind_vec(upperRightCorner_vec))   = 1;

    %% v-velocity Component Boundary Condition
    vRHS_cvec(ind_vec(lowerWall_vec))          = 0;
    vRHS_cvec(ind_vec(upperWallNoCorners_vec)) = 0;
    vRHS_cvec(ind_vec(rightWallNoCorner_vec))  = 0;
    vRHS_cvec(ind_vec(leftWallNoCorner_vec))   = 0;
    vRHS_cvec(ind_vec(upperLeftCorner_vec))    = 0;
    vRHS_cvec(ind_vec(upperRightCorner_vec))   = 0;

    vMLV_cvec(ind_vec(lowerWall_vec))          = 1;
    vMLV_cvec(ind_vec(upperWallNoCorners_vec)) = 1;
    vMLV_cvec(ind_vec(rightWallNoCorner_vec))  = 1;
    vMLV_cvec(ind_vec(leftWallNoCorner_vec))   = 1;
    vMLV_cvec(ind_vec(upperLeftCorner_vec))    = 1;
    vMLV_cvec(ind_vec(upperRightCorner_vec))   = 1;

    %% u Solution
    uSolution_cvec = uRHS_cvec./uMLV_cvec;

    %% v Solution
    vSolution_cvec = vRHS_cvec./vMLV_cvec;

    %% Reassigning initial Conditions
    un_cvec  = uSolution_cvec;
    vn_cvec  = vSolution_cvec;

    if mod(itr, vtkOutputInterval) == 0
        writeVTKatInterval(outputFolder,itr,xCoord_vec,yCoord_vec,un_cvec,vn_cvec,connectivityMatrix_mat,pressureSolution_cvec,time)
    end

    %% Stepping
    time = time + timeStep;
    itr  = itr + 1;
    disp(itr)
end
toc
writeFinalVTK(xCoord_vec,yCoord_vec,un_cvec,vn_cvec,connectivityMatrix_mat,pressureSolution_cvec,time,outputFolder)

%% Plotting Code
% variation of x-component of velocity
pakdel_mat   = readmatrix('x_velocityVariationPakdel100Re.csv', 'FileType', 'text', 'Range', 'A1:B94');
xpakdel_vec  = pakdel_mat(:,1);
ypakdel_vec  = pakdel_mat(:,2);

% Extract data at x = 0.5 and sort by y-coordinate
xEqualHalf_indices = ind_vec(xEqualHalf);
y_coords_at_xHalf  = yCoord_vec(xEqualHalf_indices);
u_vals_at_xHalf    = un_cvec(xEqualHalf_indices);

% Sort the data by y-coordinate
[sorted_y, sort_order] = sort(y_coords_at_xHalf);
sorted_u               = u_vals_at_xHalf(sort_order);

figure;
plot(sorted_u, sorted_y, 'b-', 'LineWidth', 1.5)
hold on
plot(xpakdel_vec, ypakdel_vec, 'o', 'LineWidth', 1.5)
title('variation of x-velocity component along y at x = 0.5','interpreter','latex')
xlabel('u','Interpreter','latex');
ylabel('Y','Interpreter','latex');
legend('Current Solution', 'Reference Data', 'Interpreter', 'latex')
grid on

% Change Coordinates for Re and CFL
X = reshape(xCenter_cvec, numElementsX,numElementsY);
Y = reshape(yCenter_cvec, numElementsX,numElementsY);
% Re Resutls Plotting
figure;
ReLocal_mat = reshape(ReLocal_cvec, numElementsX,numElementsY);
contourf(X, Y, ReLocal_mat, 100, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('Re Contour','Interpreter','latex');
grid on;

% CFL Resutls Plotting
figure;
CFL_mat = reshape(CFL_cvec, numElementsX,numElementsY);
contourf(X, Y, CFL_mat, 100, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('CFL Contour','Interpreter','latex');
grid on;

% % % Mesh plotting
% figure; hold on;
% X = xCoord_vec(connectivityMatrixGeo_mat);
% Y = yCoord_vec(connectivityMatrixGeo_mat);
%
% Xn = [X, nan(size(X,1),1)];
% Yn = [Y, nan(size(Y,1),1)];
%
% Xv = Xn';
% Xv = Xv(:);
% Yv = Yn';
% Yv = Yv(:);
%
% patch(Xv, Yv, 'w', 'EdgeColor', 'k');
% axis equal;
