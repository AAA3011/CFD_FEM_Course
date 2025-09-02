function [connectivityMatrix_mat,xCoord_vec,yCoord_vec] = meshContraction(lengthX_vec,lengthY_vec,numElementsX_vec,numElementsY_vec,xMin_vec,yMin_vec,totNumNodes_vec)

    numMeshes = length(lengthX_vec);
    
    connectivityMatrix_mat = [];
    xCoord_vec = [];
    yCoord_vec = [];
    
    nodeOffset = 0;
    for i = 1:numMeshes
    
        [conn,xCoord,yCoord] = structuredMesh(numElementsX_vec(i),numElementsY_vec(i),lengthX_vec(i),lengthY_vec(i),xMin_vec(i),yMin_vec(i));
        conn = conn + nodeOffset;
    
        connectivityMatrix_mat = [connectivityMatrix_mat; conn];
        xCoord_vec             = [xCoord_vec, xCoord];
        yCoord_vec             = [yCoord_vec, yCoord];
    
        nodeOffset = nodeOffset + totNumNodes_vec(i);
    end
    
    figure;
    nodes = [connectivityMatrix_mat connectivityMatrix_mat(:,1)];
    plot(xCoord_vec(nodes)', yCoord_vec(nodes)', 'k')
    title('Mesh of Connected Rectangular Domains','Interpreter','latex')
    xlabel('X','Interpreter','latex')
    ylabel('Y','Interpreter','latex')

end
