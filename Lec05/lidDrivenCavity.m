clc;clearvars;close all;

%% Connectivity Matrix and X and y Coordinates
d = 1;
lengthX = d;
lengthY = d;
numElementsX = 30;
numElementsY = 30;
xMin = 0;
yMin = 0;
[connectivityMatrix_mat,xCoord_vec,yCoord_vec] = structuredMesh(numElementsX,numElementsY,lengthX,lengthY,xMin,yMin);

xNonCoord_vec = xCoord_vec/d;
yNonCoord_vec = yCoord_vec/d;

%% Local, Global, velocity, (u,v), RHS, uRHS , vRHS Matrices

const = 0;

totNumNodes = length(xNonCoord_vec);
RHS_cvec    = zeros(totNumNodes,1);

globalMatrix_mat    = zeros(totNumNodes);
velGlobalMatrix_mat = zeros(totNumNodes);

uRHSMultAsse_mat    = zeros(totNumNodes);
vRHSMultAsse_mat    = zeros(totNumNodes);

for elementNum = 1:size(connectivityMatrix_mat,1)

    elemenetNodes_vec = connectivityMatrix_mat(elementNum,:);
    xNodesVals_vec    = xNonCoord_vec(elemenetNodes_vec);
    yNodesVals_vec    = yNonCoord_vec(elemenetNodes_vec);

    [epsiK_mat,velK_mat,uRHSMult_mat,vRHSMult_mat] = computeLocalStiffMatrix(xNodesVals_vec,yNodesVals_vec,const);

    globalMatrix_mat(elemenetNodes_vec,elemenetNodes_vec)    = globalMatrix_mat(elemenetNodes_vec,elemenetNodes_vec)    + epsiK_mat;
    velGlobalMatrix_mat(elemenetNodes_vec,elemenetNodes_vec) = velGlobalMatrix_mat(elemenetNodes_vec,elemenetNodes_vec) + velK_mat;

    uRHSMultAsse_mat(elemenetNodes_vec,elemenetNodes_vec)    = uRHSMultAsse_mat(elemenetNodes_vec,elemenetNodes_vec)    + uRHSMult_mat;
    vRHSMultAsse_mat(elemenetNodes_vec,elemenetNodes_vec)    = vRHSMultAsse_mat(elemenetNodes_vec,elemenetNodes_vec)    + vRHSMult_mat;

end

globalMatrixBeforeBC_mat      = globalMatrix_mat;
velGlobalMatrixBeforeiBC_mat  = velGlobalMatrix_mat;
velGlobalMatrixBeforeiBC2_mat = velGlobalMatrix_mat;
%% Boundary Condition

ind_vec  = 1:totNumNodes;
middelLower_vec         = (xNonCoord_vec == lengthX/2/d) & (yNonCoord_vec == 0);
upperWithoutCorners_vec = (yNonCoord_vec == lengthY/d) & (xNonCoord_vec > 0)  & (xNonCoord_vec < lengthX/d);
corners_vec             = (yNonCoord_vec == lengthY/d) & ((xNonCoord_vec == 0) | (xNonCoord_vec == lengthX/d));
everywhereElse_vec      = ~middelLower_vec & ~upperWithoutCorners_vec & ~corners_vec ;

RHS_cvec(:,1)    = 0;

globalMatrix_mat(:,:) = 0;
globalMatrix_mat(:,:) = eye(totNumNodes);

%% Stream Function Solution
streamFunctionSolution_cvec = globalMatrix_mat \ RHS_cvec;

%% RHS and solution of X-velocity Component
uRHS_cvec = uRHSMultAsse_mat * streamFunctionSolution_cvec;

uRHS_cvec(ind_vec(upperWithoutCorners_vec)) = 1;
uRHS_cvec(ind_vec(everywhereElse_vec))      = 0;
uRHS_cvec(ind_vec(corners_vec))             = 0.5;

velGlobalMatrix_mat(ind_vec(upperWithoutCorners_vec),:) = 0;
velGlobalMatrix_mat(ind_vec(upperWithoutCorners_vec),ind_vec(upperWithoutCorners_vec)) = eye(length(ind_vec(upperWithoutCorners_vec)));

velGlobalMatrix_mat(ind_vec(corners_vec),:) = 0;
velGlobalMatrix_mat(ind_vec(corners_vec),ind_vec(corners_vec)) = eye(length(ind_vec(corners_vec)));

velGlobalMatrix_mat(ind_vec(everywhereElse_vec),:) = 0;
velGlobalMatrix_mat(ind_vec(everywhereElse_vec),ind_vec(everywhereElse_vec)) = eye(length(ind_vec(everywhereElse_vec)));

uSolution_cvec = velGlobalMatrix_mat \ uRHS_cvec;

%% RHS and solution of Y-velocity Component
vRHS_cvec      = -1 * vRHSMultAsse_mat * streamFunctionSolution_cvec;

vRHS_cvec(:)    = 0;

velGlobalMatrixBeforeiBC_mat(:,:) = 0;
velGlobalMatrixBeforeiBC_mat(:,:) = eye(totNumNodes);

vSolution_cvec = velGlobalMatrixBeforeiBC_mat \ vRHS_cvec;

%% Divergence RHS and Solution
divRHS_cvec      = vRHSMultAsse_mat * uSolution_cvec + uRHSMultAsse_mat * vSolution_cvec ;
divSolution_cvec = velGlobalMatrixBeforeiBC2_mat \ divRHS_cvec;

%% Pressure Calculation
epslon = 0.05;
pressureRHS_cvec      = (-1/epslon)*(vRHSMultAsse_mat * uSolution_cvec + uRHSMultAsse_mat * vSolution_cvec) ;

pressureRHS_cvec(ind_vec(middelLower_vec))  = 0;

globalMatrixBeforeBC_mat(ind_vec(middelLower_vec),:) = 0;
globalMatrixBeforeBC_mat(ind_vec(middelLower_vec),ind_vec(middelLower_vec)) = eye(length(ind_vec(middelLower_vec)));

pressureSolution_cvec = globalMatrixBeforeBC_mat \ pressureRHS_cvec;

%% Plotting Code

%grid generation for Plotting
x_min = min(xNonCoord_vec);
x_max = max(xNonCoord_vec);
y_min = min(yNonCoord_vec);
y_max = max(yNonCoord_vec);

[x_grid, y_grid] = meshgrid(linspace(x_min, x_max, 200), linspace(y_min, y_max, 200));

% Psi Resutls Plotting
F_psi = scatteredInterpolant(xNonCoord_vec', yNonCoord_vec', streamFunctionSolution_cvec);
psi_grid = F_psi(x_grid, y_grid);

figure;
contourf(x_grid, y_grid, psi_grid, 300, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('Stream Function (\\Psi) Contour','Interpreter','latex');
grid on

figure;
contour(x_grid, y_grid, psi_grid, 300, 'LineWidth', 1.0);
% axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('Streamlines','Interpreter','latex');
grid on

% X-Velocity Component Resutls Plotting
F_u = scatteredInterpolant(xNonCoord_vec', yNonCoord_vec', uSolution_cvec);
u_grid = F_u(x_grid, y_grid);

figure;
contourf(x_grid, y_grid, u_grid, 300, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('X-Velocity Component Contour','Interpreter','latex');
grid on;

% Y-Velocity Component Resutls Plotting
F_v = scatteredInterpolant(xNonCoord_vec', yNonCoord_vec', vSolution_cvec);
v_grid = F_v(x_grid, y_grid);

figure;
contourf(x_grid, y_grid, v_grid, 300, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('Y-Velocity Component Contour','Interpreter','latex');
grid on;

% Velocity Resutls Plotting
velSolution_cvec = sqrt(uSolution_cvec.^2 + vSolution_cvec.^2);
F_vel = scatteredInterpolant(xNonCoord_vec', yNonCoord_vec', velSolution_cvec);
vel_grid = F_vel(x_grid, y_grid);

figure;
contourf(x_grid, y_grid, vel_grid, 300, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('Velocity Contour','Interpreter','latex');
grid on;

% Divergence Resutls Plotting
F_div = scatteredInterpolant(xNonCoord_vec', yNonCoord_vec', divSolution_cvec);
vel_grid = F_div(x_grid, y_grid);

figure;
contourf(x_grid, y_grid, vel_grid, 300, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('Divergenc Contour','Interpreter','latex');
grid on;

% pressure Resutls Plotting
F_p = scatteredInterpolant(xNonCoord_vec', yNonCoord_vec', pressureSolution_cvec);
p_grid = F_p(x_grid, y_grid);

figure;
contourf(x_grid, y_grid, p_grid, 300, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('Pressure Contour','Interpreter','latex');
grid on;
