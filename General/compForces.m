% compForces - Comp Forces.
%
% FILE: compForces.m
% DESCRIPTION:
% Computes drag and lift forces on cylinder by integrating stress tensor
% (pressure + viscous) over cylinder surface. Loops over boundary edges,
% averages stress at nodes, computes normal, and sums force contributions.
%
% Inputs:
%   uxSolution_col:Nodal ux-velocity solution values for the cylinder flow field.
%   uySolution_col:Nodal uy-velocity solution values for the cylinder flow field.
%   vxSolution_col:Nodal vx-velocity solution values for the cylinder flow field.
%   vySolution_col:Nodal v-velocity solution values for the cylinder flow field.
%   pressureSolution_col: Nodal pressure column vector (N x 1)
%   x_col: Nodal x-coordinates of the mesh.
%   y_col: Nodal y-coordinates of the mesh.
%   cylinderWall_vec : Node indices on the cylinder wall boundary..
%   Re (variable): Reynolds number (scalar)
% Outputs:
%   Fx : Computed total drag force in the x-direction.
%   Fy : Computed total lift force in the y-direction.
function  [Fx,Fy] = compForces(uxSolution_col,uySolution_col,vxSolution_col,vySolution_col,pressureSolution_col,x_col,y_col,cylinderWall_vec,Re)
    Fx_col  = zeros(length(cylinderWall_vec)-1,1);
    Fy_col  = zeros(length(cylinderWall_vec)-1,1);

    for i = 1:length(cylinderWall_vec)-1
        node1     = cylinderWall_vec(i);
        node2     = cylinderWall_vec(i+1);

        sigma11_1 = -pressureSolution_col(node1) + (2./Re).*uxSolution_col(node1);
        sigma12_1 = (1./Re).*(uySolution_col(node1) + vxSolution_col(node1));
        sigma21_1 = sigma12_1;
        sigma22_1 = -pressureSolution_col(node1) + (2./Re).*vySolution_col(node1);

        sigma11_2 = -pressureSolution_col(node2) + (2./Re).*uxSolution_col(node2);
        sigma12_2 = (1./Re).*(uySolution_col(node2) + vxSolution_col(node2));
        sigma21_2 = sigma12_2;
        sigma22_2 = -pressureSolution_col(node2) + (2./Re).*vySolution_col(node2);

        sigma11_avg = (sigma11_1+sigma11_2)/2;
        sigma12_avg = (sigma12_1+sigma12_2)/2;
        sigma21_avg = (sigma21_1+sigma21_2)/2;
        sigma22_avg = (sigma22_1+sigma22_2)/2;

        x1 = x_col(node1);
        x2 = x_col(node2);
        y1 = y_col(node1);
        y2 = y_col(node2);

        Le = sqrt((x2-x1)^2 + (y2-y1)^2);

        nx = (y2-y1)/Le;
        ny = -(x2-x1)/Le;

        Fx_col(i) = (sigma11_avg * nx + sigma12_avg * ny)*Le;
        Fy_col(i) = (sigma21_avg * nx + sigma22_avg * ny)*Le;
    end
    Fx = sum(Fx_col);
    Fy = sum(Fy_col);
end
