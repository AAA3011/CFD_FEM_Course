tic
clc; close all; clearvars;

%% Connectivity Matrix
connectivityMatrix_mat = [1 2;2 3;3 4;2 5;2 6;1 10;4 7;
                          5 8;6 9;7 11;8 12;9 13;10 14;
                          7 8;8 9;9 10];

%% Global Matrix parameters
numElements = size(connectivityMatrix_mat,1);
numNodes = 14;
mu = 1e-5;                                 % Dynamic Viscosity   (N * s / m^2)
d  = 0.1;                                  % Pipe Diameter       (m)
L_vec = ones(1,numElements)*10;            % Pipe Length (m) 
L_vec(2) = 20;  L_vec(6) = 20;  L_vec(4) = 14.14;

% Hydraulic conductance for each pipe
K_vec = pi * d^4 ./ (128 * L_vec * mu);

%% --------- Vectorized Assembly (no loop) ---------
node1 = connectivityMatrix_mat(:,1);
node2 = connectivityMatrix_mat(:,2);

% Row indices for all contributions
rows = [node1; node2; node1; node2];
cols = [node1; node2; node2; node1];
rows
cols
% Values for each contribution
vals = [ K_vec(:);  K_vec(:);  -K_vec(:);  -K_vec(:) ];

% Assemble sparse matrix
globalMatrix_mat = sparse(rows, cols, vals, numNodes, numNodes);

% If you need a full matrix
% globalMatrix_mat = full(globalMatrix_mat);

%% Save for later use before BC
globalMatrixBeforeBC_mat = globalMatrix_mat;

%% Boundary Conditions
RHS_cvec = zeros(numNodes,1);

% Flow rate BCs
RHS_cvec(4) = -0.71428;
RHS_cvec(6) = -0.71428;
RHS_cvec(5) = -0.71428;
RHS_cvec(1) = 5;

% Pressure BCs
pressureBCValues_vec = [1 1 1 1];           % (bar)
pressureBCNodes_vec  = [11 12 13 14];       % Node indices

for i = 1:length(pressureBCNodes_vec)
    nodeNum = pressureBCNodes_vec(i);
    temp    = globalMatrix_mat(nodeNum,nodeNum);
    globalMatrix_mat(nodeNum,:) = 0;
    RHS_cvec(:,1) = RHS_cvec(:,1) - globalMatrix_mat(:,nodeNum) * pressureBCValues_vec(i);
    globalMatrix_mat(:,nodeNum) = 0;
    globalMatrix_mat(nodeNum,nodeNum) = temp;
    RHS_cvec(nodeNum,1) = temp * pressureBCValues_vec(i);
end

%% Solve
pressureSolution_cvec = globalMatrix_mat \ RHS_cvec;
volumeFlowRate_cvec   = globalMatrixBeforeBC_mat * pressureSolution_cvec;
toc
