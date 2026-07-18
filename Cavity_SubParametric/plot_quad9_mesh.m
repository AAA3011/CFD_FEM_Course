function plot_quad9_mesh(quad9_conn,xCoord9, yCoord9)
figure; hold on; axis equal; grid on;
    
    numElems = size(quad9_conn,1);
    
    for e = 1:numElems
        conn = quad9_conn(e,:);
        
        % Plot only the corner nodes 1-2-3-4-1 to draw element edges
        order = [1 2 3 4 1];
        xe = xCoord9(conn(order));
        ye = yCoord9(conn(order));
        
        % Draw element edges
        plot(xe, ye, 'k-'); 
        
        % Draw all nodes of the element in black
        plot(xCoord9(conn), yCoord9(conn), 'k-', 'MarkerSize', 10);
        
        conn = quad9_conn(e,order(1:end-1));  % corner nodes only
        x = xCoord9(conn);
        y = yCoord9(conn);
        
        % Polygon signed area formula
        A = 0.5 * ( x(1)*(y(2)-y(4)) + x(2)*(y(3)-y(1)) + x(3)*(y(4)-y(2)) + x(4)*(y(1)-y(3)) );
        
        if A <= 0
            fprintf('Element %d has WRONG node order (CW or inverted)\n', e);
        end
    end
end
