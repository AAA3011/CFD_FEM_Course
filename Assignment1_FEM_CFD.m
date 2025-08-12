clc;close all;clearvars;

%% Connectivity Matrix

connectivityMatrix_mat = [1 2;2 3;3 4;2 5;2 6;1 10;4 7;    % Connectivity Matrix Formulation 
                          5 8;6 9;7 11;8 12;9 13;10 14;
                          7 8;8 9;9 10];

%% Global Matrix

numElements = size(connectivityMatrix_mat,1);
numNodes = 14;
mu = 10^-5 ;                                              % Dynamic Viscosity   (N * S / m^2)
d  = 0.1;                                                 % Pip Diameter        (m)
L_vec    = ones(1,numElements)*10;                        % Pip Length Initialization (m) 
L_vec(2) = 20  ;  L_vec(6) = 20   ;   L_vec(4) = 14.14;   % Pip Length Initialization (m)
K_vec    = pi * d^4 ./ (128 * L_vec * mu);                % Calculate hydraulic conductance for each pipe
globalMatrix_mat = zeros(numNodes,numNodes);              % Global Matrix Initialization 

for i = 1:numElements

    Node1 = connectivityMatrix_mat(i,1);                  % First  node of each element
    Node2 = connectivityMatrix_mat(i,2);                  % Second node of each element
    
    % Assigning values to the global Matrix 
    globalMatrix_mat(Node1,Node1) = globalMatrix_mat(Node1,Node1) + K_vec(i);
    globalMatrix_mat(Node2,Node2) = globalMatrix_mat(Node2,Node2) + K_vec(i);
    globalMatrix_mat(Node1,Node2) = globalMatrix_mat(Node1,Node2) - K_vec(i);
    globalMatrix_mat(Node2,Node1) = globalMatrix_mat(Node2,Node1) - K_vec(i);

end

globalMatrixBeforeBC_mat =  globalMatrix_mat;            % Global Marix before modifing it to fit the boundary condtions 

%% Boundary Conditions

RHS_cvec    = zeros(numNodes,1);                         % Initialize the right-hand side vector 

RHS_cvec(4) = -0.71428  ;  RHS_cvec(6) = -0.71428;       % Assigning volume flow rate boundary conditions (m^3 / S)
RHS_cvec(5) = -0.71428  ;  RHS_cvec(1) = 5;              % Assigning volume flow rate boundary conditions (m^3 / S)
 
pressureBCValues_vec = [1 1 1 1];                        % Assigning pressure boundary condition values  (bar)
pressureBCNodes_vec  = [11 12 13 14];                    % Assigning pressure boundary condition nodes

for i = 1:length(pressureBCNodes_vec)

    nodeNum = pressureBCNodes_vec(i);                                                           % get the pressure boundary condition node Num
    temp    = globalMatrix_mat(nodeNum,nodeNum);                                                % Store the diagonal element before modifying the matrix
    globalMatrix_mat(nodeNum,:) = 0;                                                            % Set all elements in the row to zero
    RHS_cvec(:,1) = RHS_cvec(:,1) - globalMatrix_mat(:,nodeNum) * pressureBCValues_vec(i);      % Update RHS vector 
    globalMatrix_mat(:,nodeNum) = 0;                                                            % Set all elements in the column to zero
    globalMatrix_mat(nodeNum,nodeNum) = temp ;                                                  % reassigning the value of the diagonal element in the global matrix
    RHS_cvec(nodeNum,1) = temp * pressureBCValues_vec(i);

end

pressureSolution_cvec = globalMatrix_mat \ RHS_cvec
volumeFlowRate_cvec   = globalMatrixBeforeBC_mat * pressureSolution_cvec 




