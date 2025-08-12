clc;clearvars;close all;


%% Connectivity Matrix

numberOfElements = 3;
connectivity_mat = zeros(numberOfElements,2);

for  i = 1:numberOfElements
    
    connectivity_mat(i,1) = i;         % Node #1
    connectivity_mat(i,2) = i+1;       % Node #2

end

%% Global Matrix

K_vec = [0.8 0.26 0.8];                                % Thermal Conductivity Vector (W/m·K) 
A_vec = [1 1 1];                                       % Cross-sectional area (m²)
h_vec = [0.01 0.02 0.01];                              % Wall thikness (m)
M_vec = zeros(1,numberOfElements);                     % M = K * A / h 
numberOfNodes = numberOfElements + 1;                  % Total Number Of Elements
globalMatrix_mat = zeros(numberOfNodes,numberOfNodes); % Globla Matrix Initialization

for i = 1:numberOfElements
   
    M_vec(i) = K_vec(i) * A_vec(i) / h_vec(i) ;
    N_1      = connectivity_mat(i,1);
    N_2      = connectivity_mat(i,2);
    globalMatrix_mat(N_1,N_1) = globalMatrix_mat(N_1,N_1) + M_vec(i) ; 
    globalMatrix_mat(N_2,N_2) = globalMatrix_mat(N_2,N_2) + M_vec(i) ; 
    globalMatrix_mat(N_1,N_2) = globalMatrix_mat(N_1,N_2) - M_vec(i) ; 
    globalMatrix_mat(N_2,N_1) = globalMatrix_mat(N_2,N_1) - M_vec(i) ; 

end

globalMatrixBeforeBC_mat =  globalMatrix_mat;

 %% Boundary Conditions

boundaryConditionsValues_vec = [233,293];              % Boundary Conditions Values4
boundaryConditionsNodes_vec  = [1,4];                  % Boundary Conditions Nodes
rightHandSide_cvec     = zeros(numberOfNodes,1); 

for  i = 1:length(boundaryConditionsNodes_vec)
    
    nodeNumber = boundaryConditionsNodes_vec(i);
    temp = globalMatrix_mat(nodeNumber,nodeNumber);   % temporary value holding 
    rightHandSide_cvec(:,1) = rightHandSide_cvec(:,1) + abs(globalMatrix_mat(:,nodeNumber)) * boundaryConditionsValues_vec(i) ;
    globalMatrix_mat(:,nodeNumber)          = 0;
    globalMatrix_mat(nodeNumber,:)          = 0;
    globalMatrix_mat(nodeNumber,nodeNumber) = temp;

end

temperatureSolution_cvec = globalMatrix_mat \ rightHandSide_cvec
heatFlux_cvec            = globalMatrixBeforeBC_mat * temperatureSolution_cvec