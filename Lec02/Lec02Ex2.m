clc;clearvars;close all;

%% Connectivity Matrix

numElements            = 100;
numNodes               = numElements+1;
connectivityMatrix_mat = zeros(numElements,2);

for i = 1:numElements
    
    connectivityMatrix_mat(i,1) = i;
    connectivityMatrix_mat(i,2) = i+1;

end

totalLength = 3;
elementLength = totalLength / numElements;
X_vec      = zeros(1,numNodes);
area_vec   = zeros(1,numNodes);
area_vec(1) = 1;

for i = 1:numElements

    X_vec(1,i+1) = X_vec(1,i) + elementLength ;
    area_vec(1,i+1) = 1 - 0.2 * X_vec(1,i);
end

%% Global Matrix

K_vec = zeros(1,numElements);
globalMatrix_mat = zeros(numNodes,numNodes);

for i = 1:numElements

    K_vec(i) = (X_vec(i+1) - X_vec(i) - 0.1 * X_vec(i+1)^2 - 0.1 * X_vec(i) ^2) / elementLength;

    globalMatrix_mat(i+1,i+1) = globalMatrix_mat(i+1,i+1) + K_vec(i);
    globalMatrix_mat(i,i)     = globalMatrix_mat(i,i) + K_vec(i);
    globalMatrix_mat(i+1,i)   = globalMatrix_mat(i+1,i) - K_vec(i);
    globalMatrix_mat(i,i+1)   = globalMatrix_mat(i,i+1) - K_vec(i);

end

globalMatrixBeforeBC_mat = globalMatrix_mat;

%% Boundary Conditions

phiDashXeq1 = 0.2;
phiX0     = 0;

RHS_cvec = zeros(numNodes,1);

RHS_cvec(end) = RHS_cvec(1) + area_vec(end) * phiDashXeq1;
RHS_cvec(1)   = phiX0;

globalMatrix_mat(1,:) = 0;
globalMatrix_mat(1,1) = 1;

phiConstants_cvec = globalMatrix_mat \ RHS_cvec;
RHSWithFirstRow   = globalMatrixBeforeBC_mat * phiConstants_cvec;

phiDashXeq0       = RHSWithFirstRow(1) / -1 * area_vec(1);

u_cvec = zeros(numNodes,1);
u_cvec(end) = phiDashXeq1 ; 

for i = 1 : numElements

    u_cvec(i) = (phiConstants_cvec(i+1) - phiConstants_cvec(i)) /  elementLength;

end

% Velocity Using FEM
uGlobalMatrix_mat = zeros(numNodes,numNodes);
uRHS_cvec         = zeros(numNodes,1);

for i = 1:numElements

    uGlobalMatrix_mat(i+1,i+1) = uGlobalMatrix_mat(i+1,i+1) + elementLength / 3;
    uGlobalMatrix_mat(i,i)     = uGlobalMatrix_mat(i,i)     + elementLength / 3;
    uGlobalMatrix_mat(i+1,i)   = uGlobalMatrix_mat(i+1,i)   + elementLength / 6;
    uGlobalMatrix_mat(i,i+1)   = uGlobalMatrix_mat(i,i+1)   + elementLength / 6;

    uRHS_cvec(i)   = uRHS_cvec(i) + (phiConstants_cvec(i+1) - phiConstants_cvec(i)) / 2 ;
    uRHS_cvec(i+1) = uRHS_cvec(i+1) + (phiConstants_cvec(i+1) - phiConstants_cvec(i)) / 2 ;

end

velocityFEM_cvec = uGlobalMatrix_mat \ uRHS_cvec;

figure(1)
plot(X_vec,area_vec,'k','LineWidth', 4)
title('Convergent Nozzle shape','interpreter','latex')
ylabel('Area $$(m^2)$$','interpreter','latex')
ylim([0 max(area_vec)* 1.3])
grid on

figure(2)
plot(X_vec,phiConstants_cvec)
title('Potential Function','interpreter','latex')
ylabel('$$\phi(x)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on


figure(3)
plot(X_vec,velocityFEM_cvec)
title('Velocity Variation Using (SG FEM)','interpreter','latex')
ylabel('Velocity $$(m/s)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on

