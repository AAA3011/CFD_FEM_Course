clc;clearvars;close all;

%% Connectivity Matrix and X and y Coordinates
Case = 'b';
d = 2;
numElements1X = 30;
numElements1Y = 30;

tolerance = 1e-10;
lengthX_vec        = [2*d d];
lengthY_vec        = [d d*0.1];
numElementsX_vec   = [numElements1X floor((lengthX_vec(2)/lengthX_vec(1))*numElements1X)];
numElementsY_vec   = [numElements1Y floor((lengthY_vec(2)/lengthY_vec(1))*numElements1Y)];
xMin_vec           = [0 lengthX_vec(1)+tolerance];
yMin_vec           = [0 0];
totNumNodes_vec    = (numElementsX_vec+1) .* (numElementsY_vec+1);

[connectivityMatrix_mat,xCoord_vec,yCoord_vec] = meshContraction(lengthX_vec,lengthY_vec,numElementsX_vec,numElementsY_vec,xMin_vec,yMin_vec,totNumNodes_vec);

xNonCoord_vec = xCoord_vec/d;
yNonCoord_vec = yCoord_vec/d;

%% Local, Global, velocity, (u,v), RHS, uRHS , vRHS Matrices

if Case == 'a'
    const = 0;
elseif Case == 'b'
    const = -(pi^2)/4;
end

totNumNodes = sum(totNumNodes_vec);
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

%% Boundary Condition

ind_vec      = 1:totNumNodes;
symmetry_vec = yNonCoord_vec == 0;
inlet_vec    = xNonCoord_vec == 0;
upper_vec    = (yNonCoord_vec == 1) | (yNonCoord_vec == 1/10 & xNonCoord_vec>= 2) | ((xNonCoord_vec == 2 & (yNonCoord_vec >= 1/10 & yNonCoord_vec <= 1)));

if Case == 'a'
    RHS_cvec(ind_vec(inlet_vec))  =  yNonCoord_vec(ind_vec(inlet_vec));
elseif Case == 'b'
    RHS_cvec(ind_vec(inlet_vec))  = sin((pi*yNonCoord_vec(ind_vec(inlet_vec)))/2);
end

RHS_cvec(ind_vec(symmetry_vec)) = 0;
RHS_cvec(ind_vec(upper_vec))    = 1;

globalMatrix_mat(ind_vec(inlet_vec),:) = 0;
globalMatrix_mat(ind_vec(inlet_vec),ind_vec(inlet_vec)) = eye(length(ind_vec(inlet_vec)));

globalMatrix_mat(ind_vec(symmetry_vec),:) = 0;
globalMatrix_mat(ind_vec(symmetry_vec),ind_vec(symmetry_vec)) = eye(length(ind_vec(symmetry_vec)));

globalMatrix_mat(ind_vec(upper_vec),:) = 0;
globalMatrix_mat(ind_vec(upper_vec),ind_vec(upper_vec)) = eye(length(ind_vec(upper_vec)));

%% Stream Function Solution
streamFunctionSolution_cvec = globalMatrix_mat \ RHS_cvec;

%% RHS and solution of X-velocity Component
uRHS_cvec      = uRHSMultAsse_mat * streamFunctionSolution_cvec;
uSolution_cvec = velGlobalMatrix_mat \ uRHS_cvec;

%% RHS and solution of Y-velocity Component
vRHS_cvec      = -vRHSMultAsse_mat * streamFunctionSolution_cvec;
vSolution_cvec = velGlobalMatrix_mat \ vRHS_cvec;

%% Divergence RHS and Solution
divRHS_cvec      = vRHSMultAsse_mat * uSolution_cvec + uRHSMultAsse_mat * vSolution_cvec ;
divSolution_cvec = velGlobalMatrix_mat \ divRHS_cvec;

%% Plotting Code

%grid generation for Plotting
x_min = min(xNonCoord_vec);
x_max = max(xNonCoord_vec);
y_min = min(yNonCoord_vec);
y_max = max(yNonCoord_vec);

[x_grid, y_grid] = meshgrid(linspace(x_min, x_max, 500), linspace(y_min, y_max, 500));
mask = (x_grid <= 2) | ((x_grid > 2) & (y_grid <= 0.1));

% Psi Resutls Plotting
F_psi = scatteredInterpolant(xNonCoord_vec', yNonCoord_vec', streamFunctionSolution_cvec);
psi_grid = F_psi(x_grid, y_grid);
psi_grid(~mask) = NaN;

figure;
contourf(x_grid, y_grid, psi_grid, 100, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('Stream Function $$(\Psi)$$ Contour','Interpreter','latex');
grid on

figure;
contour(x_grid, y_grid, psi_grid, 100, 'LineWidth', 1.0);
% axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('Streamlines','Interpreter','latex');
grid on

% X-Velocity Component Resutls Plotting
F_u = scatteredInterpolant(xNonCoord_vec', yNonCoord_vec', uSolution_cvec);
u_grid = F_u(x_grid, y_grid);
u_grid(~mask) = NaN;

figure;
contourf(x_grid, y_grid, u_grid, 100, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('X-Velocity Component Contour','Interpreter','latex');
grid on;

% Y-Velocity Component Resutls Plotting
F_v = scatteredInterpolant(xNonCoord_vec', yNonCoord_vec', vSolution_cvec);
v_grid = F_v(x_grid, y_grid);
v_grid(~mask) = NaN;

figure;
contourf(x_grid, y_grid, v_grid, 100, 'LineColor', 'none');
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
vel_grid(~mask) = NaN;

figure;
contourf(x_grid, y_grid, vel_grid, 100, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('Velocity Contour','Interpreter','latex');
grid on;

% Divergence Resutls Plotting
F_div = scatteredInterpolant(xNonCoord_vec', yNonCoord_vec', divSolution_cvec);
vel_grid = F_div(x_grid, y_grid);
vel_grid(~mask) = NaN;

figure;
contourf(x_grid, y_grid, vel_grid, 100, 'LineColor', 'none');
colorbar;
axis equal;
xlabel('x/d','Interpreter','latex');
ylabel('y/d','Interpreter','latex');
title('Divergenc Contour','Interpreter','latex');
grid on;


