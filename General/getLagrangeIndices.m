% getLagrangeIndices - Return Lagrange Indices.
%
% FILE: getLagrangeIndices.m
% DESCRIPTION:
% Compute index maps and reference coordinates for tensor-product
% Lagrange nodes on the reference interval/square with CCW ordering: corners,
% edges, then interior rings.
%
% Inputs:
%   order_1D (variable): Polynomial order of the 1D Lagrange basis (number of node intervals).
% Outputs:
%   ind_i : Row index array for tensor-product Lagrange shape functions.
%   ind_j : Column index array for tensor-product Lagrange shape functions.
%   all_nodes : Node index vector or coordinates
function [ind_i, ind_j, all_nodes] = getLagrangeIndices(order_1D)
% Same CCW ordering as before: corners -> edges -> interior rings (outer->inner)
% Minimal rewrite: rank-based traversal, no assigned-tracking needed.

    K = order_1D + 1;
    xi = linspace(-1,1,K);
    xi = xi([1, K, 2:K-1]);
    [sv, perm] = sort(xi);          % sv(r) = value at rank r; perm(r) = original index at rank r

    n = K^2;
    ind_i = zeros(1,n); ind_j = ind_i; all_nodes = zeros(n,2);
    idx = 1;

    function add(p,q)
        ind_i(idx)=perm(p); ind_j(idx)=perm(q);
        all_nodes(idx,:)=[sv(p) sv(q)];
        idx = idx+1;
    end

    lo = 1; hi = K;

    % corners CCW
    add(lo,lo); add(hi,lo); add(hi,hi); add(lo,hi);

    % outer edges CCW
    if K > 2
        I = lo+1:hi-1;
        for i=I,        add(i,lo);  end
        for j=I,        add(hi,j);  end
        for i=I(end:-1:1), add(i,hi); end
        for j=I(end:-1:1), add(lo,j); end
    end

    % interior rings CCW, outer->inner
    lo = 2; hi = K-1;
    while lo <= hi
        if lo == hi
            add(lo,lo);
        else
            for i=lo:hi,         add(i,lo);  end
            for j=lo+1:hi-1,     add(hi,j);  end
            for i=hi:-1:lo,      add(i,hi);  end
            for j=hi-1:-1:lo+1,  add(lo,j);  end
        end
        lo = lo+1; hi = hi-1;
    end

    if idx-1 ~= n
        error('getLagrangeIndices: Expected %d nodes, assigned %d',n,idx-1);
    end
end
