clc;clearvars;close all;

%% Connectivity Matrix
numElements = 100;
numNodes    = numElements+1;
connectivityMatrix_mat = zeros(numElements,2);

connectivityMatrix_mat(:,1) = 1:numElements;
connectivityMatrix_mat(:,2) = connectivityMatrix_mat(:,1) + 1;

%% Grid Generation
length = 1;
elementLength = length/numElements;

%% Stiffness Matrix and RHS
K_mat = zeros(numNodes,numNodes);
RHS_cvec = zeros(numNodes,1);

for i = 1:numElements

    node1 = connectivityMatrix_mat(i,1);
    node2 = connectivityMatrix_mat(i,2);
    
    K_mat(node1,node1) = K_mat(node1,node1) + 1/elementLength;
    K_mat(node2,node2) = K_mat(node2,node2) + 1/elementLength;
    K_mat(node2,node1) = K_mat(node2,node1) - 1/elementLength;
    K_mat(node1,node2) = K_mat(node1,node2) - 1/elementLength;
    
    RHS_cvec(node1)    = RHS_cvec(node1) + 5;
    RHS_cvec(node2)    = RHS_cvec(node2) - 5;
end

%% Boundary Condition
RHS_cvec(1) = 0;
K_mat(1,:) = 0;
K_mat(1,1) = 1;

%% Solution

ySolution_cvec = K_mat \ RHS_cvec;

figure;
x_vec = 0:elementLength:1;
plot(x_vec(:),ySolution_cvec)
ylim([-10 2])
grid on
title('Least Square Methode Demonstration','interpreter','latex')
xlabel('X','interpreter','latex')
xlabel('Y','interpreter','latex')












