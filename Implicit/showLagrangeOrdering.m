function showLagrangeOrdering(ind_i, ind_j, all_nodes, order_1D)
% Visualize the COMPLETE CCW ordering produced by getLagrangeIndices
%
% - Node numbers shown in generation order
% - Lines connect nodes sequentially
% - Corners / edges / interior colored differently

    n_nodes = size(all_nodes,1);

    % Identify node types
    node_type = zeros(n_nodes,1);  % 1 corner, 2 edge, 3 interior

    for k = 1:n_nodes
        if ( (ind_i(k)==min(ind_i) || ind_i(k)==max(ind_i)) && ...
             (ind_j(k)==min(ind_j) || ind_j(k)==max(ind_j)) )
            node_type(k) = 1; % corner
        elseif ( ind_i(k)==min(ind_i) || ind_i(k)==max(ind_i) || ...
                 ind_j(k)==min(ind_j) || ind_j(k)==max(ind_j) )
            node_type(k) = 2; % edge
        else
            node_type(k) = 3; % interior
        end
    end

    % --- Plot ---
    figure;
    hold on; grid on; axis equal;

    % Plot nodes by type
    plot(all_nodes(node_type==1,1), all_nodes(node_type==1,2), ...
         'ro','MarkerSize',10,'LineWidth',2);   % corners

    plot(all_nodes(node_type==2,1), all_nodes(node_type==2,2), ...
         'bo','MarkerSize',8,'LineWidth',1.5);  % edges

    plot(all_nodes(node_type==3,1), all_nodes(node_type==3,2), ...
         'go','MarkerSize',6,'LineWidth',1);    % interior

    % Draw ordering path
    for k = 1:n_nodes-1
        plot([all_nodes(k,1), all_nodes(k+1,1)], ...
             [all_nodes(k,2), all_nodes(k+1,2)], ...
             'k-','LineWidth',0.6);
    end

    % Label nodes
    for k = 1:n_nodes
        text(all_nodes(k,1), all_nodes(k,2), sprintf('%d',k), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',8,'FontWeight','bold');
    end

    xlabel('\xi');
    ylabel('\eta');
    title(sprintf('Q%d COMPLETE CCW Lagrange Node Ordering', order_1D));
    legend({'Corners','Edges','Interior'},'Location','bestoutside');

    hold off;
end
