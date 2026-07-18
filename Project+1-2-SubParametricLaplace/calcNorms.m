% calcNorms - Calc Norms.
%
% FILE: calcNorms.m
% DESCRIPTION:
% Compute discrete L2 and H1 norms (error norms) between the FEM solution
% and the analytical exact solution using Gauss quadrature over elements.
%
% Inputs:
%   uSolution_cvec (column vector): Nodal u-velocity values for the current field.
%   totNumElements (variable): Total number of elements in the mesh.
%   elementData (variable): Struct with precomputed element integration data (weights, shape functions, derivatives, physCoords)
% Outputs:
%   L2_norm : Discrete L2 error norm of the FEM solution over the domain.
%   H1_norm : Discrete H1 error norm of the FEM solution over the domain.
function [L2_norm,H1_norm] = calcNorms(uSolution_cvec,totNumElements,elementData)

    L2_integral = 0;
    H1_integral = 0;

    for eleNum = 1:totNumElements

        elementNodes_vec = elementData{eleNum}.Nodes;
        xNodesVals_vec   = elementData{eleNum}.xNodesVals_vec;
        yNodesVals_vec   = elementData{eleNum}.yNodesVals_vec;

        weights_pages      = elementData{eleNum}.weights;
        N_row_points_pages = elementData{eleNum}.N_row;
        J_det_points_elev  = elementData{eleNum}.J_det;
        dNdX_vec_pages     = elementData{eleNum}.dNdX;
        dNdY_vec_pages     = elementData{eleNum}.dNdY;

        % Physical coordinates at Gauss points: sum over geometry nodes
        xNodes_row = xNodesVals_vec';
        yNodes_row = yNodesVals_vec';
        
        % Compute x and y at each Gauss point
        x_gp_pages = sum(xNodes_row .* N_row_points_pages, 2);
        y_gp_pages = sum(yNodes_row .* N_row_points_pages, 2);
        
        % Evaluate exact solution at Gauss points directly
        uExact_gp  = (cosh(pi*y_gp_pages) - coth(pi)*sinh(pi*y_gp_pages)) .* sin(pi*x_gp_pages);
        uxExact_gp = pi * (cosh(pi*y_gp_pages) - coth(pi)*sinh(pi*y_gp_pages)) .* cos(pi*x_gp_pages);
        uyExact_gp = pi * (sinh(pi*y_gp_pages) - coth(pi)*cosh(pi*y_gp_pages)) .* sin(pi*x_gp_pages);

        % Evaluate FEM solution at Gauss points
        uSol_nodes = uSolution_cvec(elementNodes_vec)';

        uSol_gp  = sum(uSol_nodes .* N_row_points_pages, 2);
        uxSol_gp = sum(uSol_nodes .* dNdX_vec_pages, 2);
        uySol_gp = sum(uSol_nodes .* dNdY_vec_pages, 2);

        lError_gp = (uExact_gp - uSol_gp).^2;
        hError_gp = lError_gp + (uxExact_gp - uxSol_gp).^2 + (uyExact_gp - uySol_gp).^2;

        L2_integral = L2_integral + sum(lError_gp .* weights_pages .* J_det_points_elev, 3);
        H1_integral = H1_integral + sum(hError_gp .* weights_pages .* J_det_points_elev, 3);

    end

    L2_norm = sqrt(sum(L2_integral(:)));
    H1_norm = sqrt(sum(H1_integral(:)));
end
