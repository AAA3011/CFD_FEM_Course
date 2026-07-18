% computeDeltaVals - Compute delta Vals.
%
% FILE: computeDeltaVals.m
% DESCRIPTION:
% Solve for incremental updates (du,dv,dp) per time step using Newton
% iterations: assemble global system, apply BCs, and solve until converged.
%
% Inputs:
%   numElements (variable): Number of elements (scalar)
%   totNumNodes (variable): Total number of nodes (scalar)
%   un_col (column vector): Nodal u-velocity column vector (N x 1)
%   vn_col (column vector): Nodal v-velocity column vector (N x 1)
%   pn_col (column vector): Nodal pressure values for the current field.
%   elementData (variable): Struct with precomputed element integration data (weights, shape functions, derivatives, physCoords)
%   timeStep (variable): Time step size (scalar)
%   Re (variable): Reynolds number (scalar)
%   B (variable): Stabilization coefficient, typically equal to 1/Re when used for viscous scaling.
%   bc (variable): Struct of boundary condition node indices and prescribed values for u, v, and p.
% Outputs:
%   unplus_col : Nodal u-velocity values.
%   vnplus_col : Nodal v-velocity values.
%   pnplus_col : Nodal pressure values.
function [unplus_col,vnplus_col,pnplus_col] = computeDeltaVals(numElements,totNumNodes,un_col,vn_col,pn_col,elementData,timeStep,Re,B,bc)
    
    %% k values Initialization
    uk_col = un_col;
    vk_col = vn_col;
    pk_col = pn_col;
    
    %% Newton loop
    tol  = 1e-12;
    iter = 0;
    while (true)
        iter    = iter + 1;
        K_mat   = zeros(totNumNodes*3,totNumNodes*3);
        RHS_col = zeros(totNumNodes*3,1);
    
        for elementNumber = 1:numElements
    
            elementNodes_vec = elementData{elementNumber}.Nodes;
    
            ukNodes_col = uk_col(elementNodes_vec,1);
            vkNodes_col = vk_col(elementNodes_vec,1);
            pkNodes_col = pk_col(elementNodes_vec,1);
    
            unNodes_col = un_col(elementNodes_vec,1);
            vnNodes_col = vn_col(elementNodes_vec,1);
    
            [Ku_u_mat_pages,Ku_v_mat_pages,Ku_p_mat_pages,...
                Kv_u_mat_pages,Kv_v_mat_pages,Kv_p_mat_pages,...
                Kp_u_mat_pages,Kp_v_mat_pages,Kp_p_mat_pages,...
                Ru_col_pages,Rv_col_pages,Rp_col_pages]          = computeLocal(elementData{elementNumber},ukNodes_col,vkNodes_col,pkNodes_col,timeStep,unNodes_col,vnNodes_col,Re,B);
    
            % Build local system in node-wise ordering:
            % [du1 dv1 dp1 du2 dv2 dp2 ...]
            nLoc = numel(elementNodes_vec);
            idxU = 1:3:(3*nLoc);
            idxV = 2:3:(3*nLoc);
            idxP = 3:3:(3*nLoc);
    
            %% Preparing Local K and RHS
            KLocal_mat = zeros(3*nLoc, 3*nLoc);
            KLocal_mat(idxU, idxU) = Ku_u_mat_pages;
            KLocal_mat(idxU, idxV) = Ku_v_mat_pages;
            KLocal_mat(idxU, idxP) = Ku_p_mat_pages;
    
            KLocal_mat(idxV, idxU) = Kv_u_mat_pages;
            KLocal_mat(idxV, idxV) = Kv_v_mat_pages;
            KLocal_mat(idxV, idxP) = Kv_p_mat_pages;
    
            KLocal_mat(idxP, idxU) = Kp_u_mat_pages;
            KLocal_mat(idxP, idxV) = Kp_v_mat_pages;
            KLocal_mat(idxP, idxP) = Kp_p_mat_pages;
    
            RHSLocal_col = zeros(3*nLoc, 1);
            RHSLocal_col(idxU) = Ru_col_pages;
            RHSLocal_col(idxV) = Rv_col_pages;
            RHSLocal_col(idxP) = Rp_col_pages;
    
            ids = reshape([3*elementNodes_vec - 2; 3*elementNodes_vec - 1; 3*elementNodes_vec], 1, []);
    
            %% Global Matrix Assembly
            K_mat(ids, ids) = K_mat(ids, ids) + KLocal_mat;
            RHS_col(ids)    = RHS_col(ids) + RHSLocal_col;
        end
    
        [K_mat, RHS_col] = applyBCOnIncrements(K_mat, RHS_col, uk_col, vk_col, pk_col, bc);
    
        K_mat = decomposition(K_mat);
        uvpDelta_col = K_mat \ RHS_col;
        nor = norm(RHS_col, 2);
    
        uDelta_col   = uvpDelta_col(1:3:end);
        vDelta_col   = uvpDelta_col(2:3:end);
        pDelta_col   = uvpDelta_col(3:3:end);
    
        %% Calc Values after one step
        uKplus_col   = uk_col + uDelta_col;
        vKplus_col   = vk_col + vDelta_col;
        pKplus_col   = pk_col + pDelta_col;
    
        %% assign previous step's values the new values
        uk_col = uKplus_col;
        vk_col = vKplus_col;
        pk_col = pKplus_col;
    
        if (nor < tol)
            break;
        end
    end
    
    unplus_col = uKplus_col;
    vnplus_col = vKplus_col;
    pnplus_col = pKplus_col;

end
