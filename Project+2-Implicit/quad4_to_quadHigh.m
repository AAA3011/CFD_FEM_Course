% quad4_to_quadHigh - Convert Quad4 mesh to higher-order quadrilateral mesh.
%
% FILE: quad4_to_quadHigh.m
% DESCRIPTION:
% Convert a 4-node quadrilateral mesh to a higher-order quadrilateral
% representation (adds edge and interior nodes for a given polynomial order).
%
% Inputs:
%   quad4_conn (variable): Connectivity matrix
%   xCoord (variable): Vector of x-coordinates.
%   yCoord (variable): Vector of y-coordinates.
%   order (variable): Polynomial order of the solution element basis.
% Outputs:
%   quadN_conn : Connectivity matrix
%   xCoordN : Vector of x-coordinates.
%   yCoordN : Vector of y-coordinates.
function [quadN_conn, xCoordN, yCoordN] = quad4_to_quadHigh(quad4_conn, xCoord, yCoord, order)

    if nargin < 4, error('Need: quad4_conn, xCoord, yCoord, order'); end
    p = double(order);
    if p < 1 || p ~= floor(p), error('order must be integer >= 1'); end

    nE   = size(quad4_conn, 1);
    nN   = numel(xCoord);
    nPE  = (p+1)^2;                          % nodes per element

    % Pre-allocate coords, copy originals
    cap      = nN + nE*((p-1)^2 + 4*(p-1));
    xCoordN  = [xCoord(:); nan(cap - nN, 1)];
    yCoordN  = [yCoord(:); nan(cap - nN, 1)];
    nextNode = nN + 1;

    quadN_conn = zeros(nE, nPE);
    edgeMap    = containers.Map('KeyType','char','ValueType','double');

    % Bilinear map: reference (xi,eta) in [-1,1]x[-1,1] -> physical (x,y)
    % Standard shape functions, corners CCW: 1(-1,-1) 2(1,-1) 3(1,1) 4(-1,1)
    bmap = @(xi,eta, x1,x2,x3,x4, y1,y2,y3,y4) deal( ...
        0.25*(1-xi).*(1-eta).*x1 + 0.25*(1+xi).*(1-eta).*x2 + 0.25*(1+xi).*(1+eta).*x3 + 0.25*(1-xi).*(1+eta).*x4, ...
        0.25*(1-xi).*(1-eta).*y1 + 0.25*(1+xi).*(1-eta).*y2 + 0.25*(1+xi).*(1+eta).*y3 + 0.25*(1-xi).*(1+eta).*y4 );

    for e = 1:nE
        c = quad4_conn(e,:);                  % 4 corner node indices (CCW)
        quadN_conn(e,1:4) = c;
        [x1,x2,x3,x4] = deal(xCoord(c(1)), xCoord(c(2)), xCoord(c(3)), xCoord(c(4)));
        [y1,y2,y3,y4] = deal(yCoord(c(1)), yCoord(c(2)), yCoord(c(3)), yCoord(c(4)));
        pos = 5;

        % ---- Edge nodes (shared between adjacent elements) ----
        % Edge param t=k/p in (0,1), rescaled to s=2t-1 in (-1,1)
        % direction: edge1(eta=-1), edge2(xi=1), edge3(eta=1,reversed), edge4(xi=-1,reversed)
        edgePairs = [1,2; 2,3; 3,4; 4,1];
        for edge = 1:4
            nA = c(edgePairs(edge,1));  nB = c(edgePairs(edge,2));
            for k = 1:(p-1)
                t = k/p;
                s = 2*t - 1;
                switch edge
                    case 1, xi=s;   eta=-1;
                    case 2, xi=1;   eta=s;
                    case 3, xi=-s;  eta=1;
                    case 4, xi=-1;  eta=-s;
                end

                % Canonical edge key: always ordered min->max node, k_key counts from min end
                k_key = k; if nA > nB, k_key = p-k; end
                key = sprintf('%d-%d-%d', min(nA,nB), max(nA,nB), k_key);

                if isKey(edgeMap, key)
                    nid = edgeMap(key);
                else
                    [xp, yp] = bmap(xi,eta, x1,x2,x3,x4, y1,y2,y3,y4);
                    nid = nextNode;
                    xCoordN(nid) = xp;  yCoordN(nid) = yp;
                    edgeMap(key) = nid;
                    nextNode = nextNode + 1;
                end
                quadN_conn(e, pos) = nid;  pos = pos + 1;
            end
        end

        % ---- Interior nodes in CCW concentric rings ----
        m = p-1;
        if m > 0
            % Create all interior nodes on an m×m grid (i=xi-index, j=eta-index)
            ids = zeros(m,m);
            for j = 1:m
                for i = 1:m
                    xi_ij  = 2*(i/p) - 1;
                    eta_ij = 2*(j/p) - 1;
                    [xp, yp] = bmap(xi_ij, eta_ij, x1,x2,x3,x4, y1,y2,y3,y4);
                    ids(i,j) = nextNode;
                    xCoordN(nextNode) = xp;  yCoordN(nextNode) = yp;
                    nextNode = nextNode + 1;
                end
            end

            % Walk rings from outside in, CCW: bottom→right→top→left
            for r = 1:ceil(m/2)
                L=r; R=m-r+1; B=r; T=m-r+1;
                if L==R && B==T                           % single center node
                    quadN_conn(e,pos) = ids(L,B);  pos = pos+1;
                else
                    for i = L:R,        quadN_conn(e,pos) = ids(i,B); pos=pos+1; end  % bottom L→R
                    for j = B+1:T,      quadN_conn(e,pos) = ids(R,j); pos=pos+1; end  % right  B→T
                    if T~=B, for i = R-1:-1:L,  quadN_conn(e,pos) = ids(i,T); pos=pos+1; end; end  % top R→L
                    if L~=R, for j = T-1:-1:B+1, quadN_conn(e,pos) = ids(L,j); pos=pos+1; end; end % left T→B
                end
            end
        end

        assert(pos-1 == nPE, 'Node count mismatch in element %d', e);
    end

    % Trim to actual node count
    xCoordN = xCoordN(1:nextNode-1);
    yCoordN = yCoordN(1:nextNode-1);

end
