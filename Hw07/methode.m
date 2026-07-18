function [TSolution_cvec] = methode(numElements,elementLength,Pe,F,Tho)

    %% Connectivity Matrix
    numNodes    = numElements+1;
    connectivityMatrix_mat = zeros(numElements,2);
    
    connectivityMatrix_mat(:,1) = 1:numElements;
    connectivityMatrix_mat(:,2) = connectivityMatrix_mat(:,1) + 1;
    
    %% Stiffness Matrix and RHS
    K_mat    = zeros(numNodes,numNodes);
    RHS_cvec = zeros(numNodes,1);

    term = (F+Tho*Pe^2)/elementLength;
    
    for i = 1:numElements
    
        node1 = connectivityMatrix_mat(i,1);
        node2 = connectivityMatrix_mat(i,2);
    
        K_mat(node1,node1) = K_mat(node1,node1) + ( term - Pe*F/2);
        K_mat(node2,node2) = K_mat(node2,node2) + ( term + Pe*F/2);
        K_mat(node2,node1) = K_mat(node2,node1) - ( term + Pe*F/2);
        K_mat(node1,node2) = K_mat(node1,node2) + (-term + Pe*F/2);
        
    end
    
    %% Boundary Condition
    RHS_cvec(1)   = 0;
    RHS_cvec(end) = 1;
    
    K_mat(1,:) = 0;
    K_mat(1,1) = 1;
    
    K_mat(end,:)   = 0;
    K_mat(end,end) = 1;
    
    %% Solution
    TSolution_cvec = K_mat \ RHS_cvec;

end