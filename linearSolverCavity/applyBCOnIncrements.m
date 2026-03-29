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
    K_global(:,row)    = 0;
    K_global(row, row) = 1;
    R_global(row) = rhsVal;

end
