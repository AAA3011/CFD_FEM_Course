% applyBCOnIncrements - Enforce Dirichlet BCs on solution increments.
%
% FILE: applyBCOnIncrements.m
% DESCRIPTION:
% Modify the global system matrix `K_global` and RHS `R_global` to
% enforce prescribed increments (Dirichlet BCs) on velocity and
% pressure degrees of freedom for an incremental solver. This sets the
% corresponding rows to zero except the diagonal and assigns the
% target increment value into the RHS.
%
% Inputs:
%   K_global : Global stiffness/linearized system matrix (sparse or dense).
%   R_global : Global right-hand-side vector corresponding to increments.
%   u_col    : Current nodal u-velocity values (N x 1 column vector).
%   v_col    : Current nodal v-velocity values (N x 1 column vector).
%   p_col    : Current nodal pressure values (N x 1 column vector).
%   bc       : Struct with boundary condition info containing fields:
%                - uNodes, uVals : node indices and prescribed u-values
%                - vNodes, vVals : node indices and prescribed v-values
%                - pNodes, pVals : node indices and prescribed p-values
%
% Outputs:
%   K_global : Modified system matrix with BC rows enforced.
%   R_global : Modified RHS vector with BC values applied.
function [K_global, R_global] = applyBCOnIncrements(K_global, R_global, u_col, v_col, p_col, bc)

    for k = 1:length(bc.uNodes)
        node = bc.uNodes(k);
        row  = 3*node - 2;
        targetDu = bc.uVals(k) - u_col(node);

        [K_global, R_global] = imposeSingleIncrementBC(K_global, R_global, row, targetDu);
    end


    for k = 1:length(bc.vNodes)
        node = bc.vNodes(k);
        row  = 3*node - 1;
        targetDv = bc.vVals(k) - v_col(node);

        [K_global, R_global] = imposeSingleIncrementBC(K_global, R_global, row, targetDv);
    end

    for k = 1:length(bc.pNodes)
        node = bc.pNodes(k);
        row  = 3*node;
        targetDp = bc.pVals(k) - p_col(node);

        [K_global, R_global] = imposeSingleIncrementBC(K_global, R_global, row, targetDp);
    end
end

function [K_global, R_global] = imposeSingleIncrementBC(K_global, R_global, row, rhsVal)

    K_global(row, :)   = 0;
    % K_global(:,row)    = 0;
    K_global(row, row) = 1;
    R_global(row) = rhsVal;

end
