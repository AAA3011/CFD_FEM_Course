clc;clearvars;close all;
tic
%% Element informations
numGaussPoints  = 10;
numGeometryGaussPoints = numGaussPoints;

order   = 1;                % Order of Solution Element
K       = (order + 1)^2;    % Solution Element number of nodes (has to be consistent with order of element specified below)

order_G = 1;                % Order of Geometric Element (Should be the same as Geometric meshing)
K_G     = (order_G + 1)^2;  % Geometric Element number of nodes (has to be consistent with order of element specified below)

%% Global Parameters
Re        = 10;
B         = 1/Re;
timeStep  = 0.0005;
totitrs   = 30;

%% VTK file realated parameters
vtkOutputInterval = 1000;
outputFolder      = 'CavityResutts100ReHigher';
makeVTKdir(outputFolder)

%% Connectivity Matrix and Grid Generation
numElementsX = 40;
numElementsY = 40;
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

%% Precompute geometric data for all elements 
numElements = numElementsX * numElementsY;
elementData = cell(numElements, 1);
for eleNum = 1:numElements
    elementNodesGeo_vec = connectivityMatrixGeo_mat(eleNum,:);
    xNodesValsGeo_vec   = xCoordGeo_vec(elementNodesGeo_vec);
    yNodesValsGeo_vec   = yCoordGeo_vec(elementNodesGeo_vec);
    [weights_pages,N_row_points_pages,J_det_points_elev,N_diff_PhysCoords_points_rows_pages,ElementPhysCoords_mat] = getGaussRelated(xNodesValsGeo_vec,yNodesValsGeo_vec,K,numGaussPoints,K_G,numGeometryGaussPoints,order,order_G);
    elementData{eleNum}.weights    = weights_pages;
    elementData{eleNum}.N_row      = N_row_points_pages;
    elementData{eleNum}.J_det      = J_det_points_elev;
    elementData{eleNum}.dNdX       = N_diff_PhysCoords_points_rows_pages(1, :, :);
    elementData{eleNum}.dNdY       = N_diff_PhysCoords_points_rows_pages(2, :, :);
    elementData{eleNum}.physCoords = ElementPhysCoords_mat;
end

%% Boundary Conditions locations
tolBC   = 1e-12;
ind_vec = 1:length(xCoord_vec);

lowerWall_vec          =  abs(yCoord_vec - 0) < tolBC;
upperWallNoCorners_vec = (abs(yCoord_vec - lengthY) < tolBC) & (abs(xCoord_vec - lengthX) >= tolBC) & (abs(xCoord_vec - 0) >= tolBC);
rightWallNoCorner_vec  = (abs(xCoord_vec - lengthX) < tolBC) & (abs(yCoord_vec - lengthY) >= tolBC);
leftWallNoCorner_vec   = (abs(xCoord_vec - 0) < tolBC)       & (abs(yCoord_vec - lengthY) >= tolBC);
upperLeftCorner_vec    = (abs(xCoord_vec - 0) < tolBC)       & (abs(yCoord_vec - lengthY) < tolBC);
upperRightCorner_vec   = (abs(xCoord_vec - lengthX) < tolBC) & (abs(yCoord_vec - lengthY) < tolBC);
pBCCondition_col       = (abs(xCoord_vec - 0.5)< tolBC) & (abs(yCoord_vec - 0)< tolBC);

%% Initial Conditions
totNumNodes = length(xCoord_high);
un_col = zeros(totNumNodes,1);
vn_col = zeros(totNumNodes,1);
pn_col = zeros(totNumNodes,1);

%% Boundary Conditions Informations
uNodesAll = [ ...
    ind_vec(lowerWall_vec), ...
    ind_vec(upperWallNoCorners_vec), ...
    ind_vec(rightWallNoCorner_vec), ...
    ind_vec(leftWallNoCorner_vec), ...
    ind_vec(upperLeftCorner_vec), ...
    ind_vec(upperRightCorner_vec)];
    
uValsAll = [ ...
    zeros(1,nnz(lowerWall_vec)), ...
    ones(1,nnz(upperWallNoCorners_vec)), ...
    zeros(1,nnz(rightWallNoCorner_vec)), ...
    zeros(1,nnz(leftWallNoCorner_vec)), ...
    0.5, 0.5];

vNodesAll = uNodesAll;
vValsAll  = zeros(size(vNodesAll));

pNodesAll = ind_vec(pBCCondition_col);
pValsAll  = zeros(size(pNodesAll));

[bc.uNodes, iu] = unique(uNodesAll, "stable");
bc.uVals = uValsAll(iu).';

[bc.vNodes, iv] = unique(vNodesAll, "stable");
bc.vVals = vValsAll(iv).';

[bc.pNodes, ip] = unique(pNodesAll, "stable");
bc.pVals = pValsAll(ip).';

%% Enforce BC on initialized fields
un_col(bc.uNodes) = bc.uVals;
vn_col(bc.vNodes) = bc.vVals;
pn_col(bc.pNodes) = bc.pVals;

%% time Loop
itr  = 1;
time = 0;
Diver_vec = [];
while(itr < totitrs)

    %% Compute Delta Values
    [Diver,unplus_col,vnplus_col,pnplus_col] = computeDeltaVals(numElements,totNumNodes,connectivityMatrix_mat,un_col,vn_col,pn_col,elementData,timeStep,Re,B,bc);
    Diver_vec(itr) = Diver;
    %% Reassigning initial Conditions
    un_col  = unplus_col;
    vn_col  = vnplus_col;
    pn_col  = pnplus_col;

    %% Writting VTK File
    if mod(itr, vtkOutputInterval) == 0
        writeVTKatInterval(outputFolder,itr,xCoord_vec,yCoord_vec,un_col,vn_col,connectivityMatrix_mat,pn_col,time)
    end

    %% Stepping
    time = time + timeStep;
    itr  = itr + 1;
    disp(itr)
    
end
toc
writeFinalVTK(xCoord_vec,yCoord_vec,un_col,vn_col,connectivityMatrix_mat,pn_col,time,outputFolder)

%% Plotting Code
% variation of x-component of velocity
pakdel_mat   = readmatrix('x_velocityVariationPakdel100Re.csv', 'FileType', 'text', 'Range', 'A1:B94');
xpakdel_vec  = pakdel_mat(:,1);
ypakdel_vec  = pakdel_mat(:,2);

% Extract data at x = 0.5 and sort by y-coordinate
xEqualHalf_indices = ind_vec(xCoord_vec == 0.5);
y_coords_at_xHalf  = yCoord_vec(xEqualHalf_indices);
u_vals_at_xHalf    = un_col(xEqualHalf_indices);

% Sort the data by y-coordinate
[sorted_y, sort_order] = sort(y_coords_at_xHalf);
sorted_u = u_vals_at_xHalf(sort_order);

figure;
plot(sorted_u, sorted_y, 'b-', 'LineWidth', 1.5)
hold on
plot(xpakdel_vec, ypakdel_vec, 'o', 'LineWidth', 1.5)
title(sprintf(['Variation of x-velocity component along y at x = 0.5. ' ...
               'Taken Time: %.3f min, Element used: Q%d, Num of Elements: %d'], ...
               toc/60, K, numElements), ...
      'Interpreter','latex')
xlabel('u','Interpreter','latex');
ylabel('Y','Interpreter','latex');
legend('Current Solution', 'Reference Data', 'Interpreter', 'latex')
grid on

%% Contour Plots: p, u, v distributions
nPlot = 121;
xPlot = linspace(min(xCoord_vec), max(xCoord_vec), nPlot);
yPlot = linspace(min(yCoord_vec), max(yCoord_vec), nPlot);
[Xplot, Yplot] = meshgrid(xPlot, yPlot);

Fu = scatteredInterpolant(xCoord_vec(:), yCoord_vec(:), un_col(:), 'natural', 'nearest');
Fv = scatteredInterpolant(xCoord_vec(:), yCoord_vec(:), vn_col(:), 'natural', 'nearest');
Fp = scatteredInterpolant(xCoord_vec(:), yCoord_vec(:), pn_col(:), 'natural', 'nearest');

Uplot = Fu(Xplot, Yplot);
Vplot = Fv(Xplot, Yplot);
Pplot = Fp(Xplot, Yplot);
NUM = 500;
figure;
contourf(Xplot, Yplot, Pplot, NUM, 'LineColor', 'none');
colorbar;
axis equal tight;
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
title('Pressure Contour Distribution','Interpreter','latex');
grid on;

figure;
contourf(Xplot, Yplot, Uplot, NUM, 'LineColor', 'none');
colorbar;
axis equal tight;
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
title('u-Velocity Contour Distribution','Interpreter','latex');
grid on;

figure;
contourf(Xplot, Yplot, Vplot, NUM, 'LineColor', 'none');
colorbar;
axis equal tight;
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
title('v-Velocity Contour Distribution','Interpreter','latex');
grid on;

%% Divervence Plotting
figure;
plot(1:totitrs-1,Diver_vec)