clc;clearvars;close all;

%% Connectivity Matrix

timeStep      = 0.001;
totalTime     = 1.1;
Time          = 0;

%% Connectivity Matrix
numElements             = 20;
numNodes                = numElements + 1;
connvectivityMatrix_mat = zeros(numElements,2);

connvectivityMatrix_mat(:,1) = 1:numElements;
connvectivityMatrix_mat(:,2) = connvectivityMatrix_mat(:,1) + 1;

%% Grid Generation
domainStart    = -1;
domainEnd      = 1;
domainLength   = abs(domainStart)+abs(domainEnd);
elementLength  = domainLength / numElements ;
x_vec          = domainStart:elementLength:domainEnd;

%% initial Condition
Tn_cvec   = 1-(x_vec.^2)';

%% Global/local Stiffness and Mass Matrix
K_mat     = zeros(numNodes,numNodes);
M_mat     = zeros(numNodes,numNodes);

for i = 1:numElements

    node1 = connvectivityMatrix_mat(i,1);
    node2 = connvectivityMatrix_mat(i,2);

    elementNodes        = [node1 node2];
    physicalCoordinates = [x_vec(node1) x_vec(node2)];

    [k_mat,m_mat] = computeLocal(physicalCoordinates);

    K_mat(elementNodes,elementNodes) = K_mat(elementNodes,elementNodes) + k_mat;
    M_mat(elementNodes,elementNodes) = M_mat(elementNodes,elementNodes) + m_mat;
end

%% Boundary Condtions and RHS_cvec

M_mat(1,:)     = 0;
M_mat(1,1)     = 1;
M_mat(end,:)   = 0;
M_mat(end,end) = 1;

while(Time  <= totalTime )

    RHS_cvec            = (M_mat - timeStep*K_mat) * Tn_cvec;
    RHS_cvec(1)         = 0;
    RHS_cvec(end)       = 0;

    Tsolution_cvec = M_mat \ RHS_cvec;
    Tn_cvec        = Tsolution_cvec;

    if any(abs(Time - [0,0.25,0.5,0.75,1]) < 1e-6)
        plot(x_vec',Tsolution_cvec,'DisplayName', [num2str(Time) ' sec'])
        title('variatio of Temperature at different time steps','Interpreter','latex')
        xlabel('X','Interpreter','latex')
        ylabel('Temperature $$(K)$$','Interpreter','latex')
        hold on
        grid on
    end
    Time = Time +timeStep;
end

grid on
legend show
