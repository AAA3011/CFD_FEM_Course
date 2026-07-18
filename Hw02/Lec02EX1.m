clc;clearvars;close all;

%% Connectivity Matrix

numElements = 3;
numNodes    = numElements + 1 ;
connectivityMatrix_mat = zeros(numElements,2);

for i =  1:numElements

    connectivityMatrix_mat(i,1) = i;
    connectivityMatrix_mat(i,2) = i+1;

end

L_vec = [0.6 0.3 0.1];
K_vec = 1./ L_vec;


%% Global Matrix 

globalMatrix_mat = zeros(numNodes);
RHS_cvec         = zeros(numNodes,1);

for i = 1 : numElements

    globalMatrix_mat(i,i)     = globalMatrix_mat(i,i)     + K_vec(i);
    globalMatrix_mat(i+1,i+1) = globalMatrix_mat(i+1,i+1) + K_vec(i);
    globalMatrix_mat(i+1,i)   = globalMatrix_mat(i+1,i)   - K_vec(i);
    globalMatrix_mat(i,i+1)   = globalMatrix_mat(i,i+1)   - K_vec(i);

    RHS_cvec(i,1)   = RHS_cvec(i,1) + L_vec(i)*0.5;
    RHS_cvec(i+1,1) = RHS_cvec(i+1,1) + L_vec(i)*0.5;

end
globalMatrixBeforeBC_mat  = globalMatrix_mat;
RHSBeforeBC_cvec          = RHS_cvec;

%% Boundary Conditions

yBCValX1   = 1;
yDashX0  = 0; 

RHS_cvec(1,1) = RHS_cvec(1,1) - yDashX0 ; 

RHS_cvec(end,end) = 0;
RHS_cvec(end,end) = yBCValX1;

globalMatrix_mat(end,:) = 0;
globalMatrix_mat(end,end) = 1;

%% Solution

nodeConsts  = globalMatrix_mat \ RHS_cvec ;
RHSWithLastRow = globalMatrixBeforeBC_mat * nodeConsts ; 

yDashX1        = RHSWithLastRow(end) - RHSBeforeBC_cvec(end);

figure(1)
plot([0 0.6 0.9 1],nodeConsts)
title('Y Function','interpreter','latex')
ylabel('$$y(x)$$','interpreter','latex')
xlabel('X','interpreter','latex')
ylim([0 max(nodeConsts)*1.2])
grid on
