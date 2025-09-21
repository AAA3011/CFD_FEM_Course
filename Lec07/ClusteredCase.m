function[TSolution_cvec,PeLocal_vec,x_vec] = ClusteredCase(numElements,length,rate,F,Tho,Pe)
    
    %% Connectivity Matrix
    numNodes    = numElements+1;
    connectivityMatrix_mat = zeros(numElements,2);
    
    connectivityMatrix_mat(:,1) = 1:numElements;
    connectivityMatrix_mat(:,2) = connectivityMatrix_mat(:,1) + 1;
    
    %% Clustered Grid Generation
    indexNodeX_vec = 1:numNodes;
    x_vec = length * (1 - rate.^(indexNodeX_vec-1)) / (1 - rate^numElements);
    
    %% Stiffness Matrix and RHS
    K_mat        = zeros(numNodes,numNodes);
    RHS_cvec     = zeros(numNodes,1);
    PeLocal_vec  = zeros(1,numNodes);
    
    for i = 1:numElements
    
        node1 = connectivityMatrix_mat(i,1);
        node2 = connectivityMatrix_mat(i,2);
    
        elementLength      = x_vec(node2)-x_vec(node1);

        if i == 1 
            PeLocal_vec(node1) = Pe*elementLength/length;
            PeLocal_vec(node2) = Pe*elementLength/length;
        else
            PeLocal_vec(node2) = Pe*elementLength/length;
            PeLocal_vec(node1) = (PeLocal_vec(node2)+PeLocal_vec(node1))/2;
        end

        if Tho == 2
            Tho = 1/sqrt((2*Pe/elementLength)^2+(4*-1/elementLength)^2);
        end
        term              = (1+Tho*Pe^2)/elementLength;
    
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
    
    %% GLS Solution
    TSolution_cvec = K_mat \ RHS_cvec;
end