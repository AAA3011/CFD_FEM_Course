clc;clearvars;close all;

%% Global Parameters
numGaussPoints = 4;
n_ElementType  = 0;
timeStep       = 0.001;
mu_vec         = [0.0001 0.001 0.0025 0.005];
totalTime      = 1.1;

%% Connectivity Matrix
numElements  = 200;
numNodes     = numElements + 1;
connvectivityMatrix_mat = zeros(numElements,2);

connvectivityMatrix_mat(:,1) = 1:numElements;
connvectivityMatrix_mat(:,2) = connvectivityMatrix_mat(:,1) + 1;

%% Grid Generation
domainStart    = -1;
domainEnd      = 1;
itr            = 0;
domainLength   = abs(domainStart)+abs(domainEnd);
elementLength  = domainLength / numElements ;
x_vec          = domainStart:elementLength:domainEnd;


for j = 1:length(mu_vec)
    %% initialization
    initialTime     = 0;
    Time            = 0;
    uShock_cvec     = zeros(numNodes,1);
    uExpantion_cvec = zeros(numNodes,1);

    %% initial Condtions
    ind = 1:length(x_vec);
    left  = x_vec >= -1 & x_vec < 0;
    Right = x_vec >=  0  & x_vec <= 1;

    un_cvec        = exp(-16 * x_vec'.^2);

    unLump_cvec    = un_cvec;
    unLumpvec_cvec = un_cvec;

    uShock_cvec(ind(left))     = 1;
    uShock_cvec(ind(Right))    = 0;

    uExpantion_cvec(ind(left))    = 0;
    uExpantion_cvec(ind(Right))   = 1;

    MExpantion_cvec  = zeros(numNodes,1);
    MShock_cvec      = zeros(numNodes,1);
    M_lump           = zeros(numNodes,numNodes);
    M_cvec           = zeros(numNodes,1);
    for i = 1:numElements

        node1 = connvectivityMatrix_mat(i,1);
        node2 = connvectivityMatrix_mat(i,2);

        elementNodes        = [node1 node2];
        physicalCoordinates = [x_vec(node1) x_vec(node2)];

        [m_lump,m_cvec,mExpantion_cvec,mShock_cvec]=compLocalTimeIndependent(physicalCoordinates,numGaussPoints,n_ElementType);

        M_cvec(elementNodes,1)            = M_cvec(elementNodes,1)            + m_cvec;
        M_lump(elementNodes,elementNodes) = M_lump(elementNodes,elementNodes) + m_lump;
        MShock_cvec(elementNodes,1)       = MShock_cvec(elementNodes,1)       + mShock_cvec;
        MExpantion_cvec(elementNodes,1)   = MExpantion_cvec(elementNodes,1)   + mExpantion_cvec;
    end
    figure;
    while(Time  <= totalTime )
        R_cvec           = zeros(numNodes,1);
        M_mat            = zeros(numNodes,numNodes);
        RShock_cvec      = zeros(numNodes,1);
        RExpantion_cvec  = zeros(numNodes,1);

        for i = 1:numElements

            node1 = connvectivityMatrix_mat(i,1);
            node2 = connvectivityMatrix_mat(i,2);

            elementNodes        = [node1 node2];
            physicalCoordinates = [x_vec(node1) x_vec(node2)];
            u_nodes             = [un_cvec(node1) un_cvec(node2)]';
            uShock_nodes        = [uShock_cvec(node1) uShock_cvec(node2)]';
            uExpantion_nodes    = [uExpantion_cvec(node1) uExpantion_cvec(node2)]';
            [r_cvec,m_mat,rShock_cvec,rExpantion_cvec] = compLocal(physicalCoordinates,u_nodes,mu_vec(j),timeStep,elementLength,uShock_nodes,uExpantion_nodes);

            R_cvec(elementNodes,1)            = R_cvec(elementNodes,1)            + r_cvec;
            M_mat(elementNodes,elementNodes)  = M_mat(elementNodes,elementNodes)  + m_mat;

            RShock_cvec(elementNodes,1)       = RShock_cvec(elementNodes,1)       + rShock_cvec;

            RExpantion_cvec(elementNodes,1)   = RExpantion_cvec(elementNodes,1)   + rExpantion_cvec;
        end

        M_lump = sum(M_lump,2);

        %% FUll Mass matrix RHS and BC
        RHS_cvec       = -timeStep * R_cvec + M_mat*un_cvec;
        RHS_cvec(1)    = exp(-16);
        RHS_cvec(end)  = exp(-16);

        M_mat(1,:)     = 0;
        M_mat(1,1)     = 1;
        M_mat(end,:)   = 0;
        M_mat(end,end) = 1;

        uSolution_cvec = M_mat \ RHS_cvec;
        un_cvec        = uSolution_cvec;

        %% Lumped Solution RHS and BC

        RHSLumped_cvec      = -timeStep * R_cvec + M_lump.*unLump_cvec;
        RHSLumped_cvec(1)   = exp(-16);
        RHSLumped_cvec(end) = exp(-16);

        M_lump(1,1)     = 1;
        M_lump(end,end) = 1;

        uSolutionLumped_mat = RHSLumped_cvec./M_lump;
        unLump_cvec          = uSolutionLumped_mat;

        %% Vector Lumping RHS and BC
        RHSLumpedcvec_cvec       = -timeStep * R_cvec + M_cvec.*unLumpvec_cvec;
        RHSLumpedcvec_cvec(1)    = exp(-16);
        RHSLumpedcvec_cvec(end)  = exp(-16);


        M_cvec(1,1)     = 1;
        M_cvec(end,end) = 1;

        uSolutionLumpedcvec_cvec = RHSLumpedcvec_cvec./M_cvec;
        unLumpvec_cvec           = uSolutionLumpedcvec_cvec;

        %% Shock Case RHS and BC
        RHSShock_cvec       = -timeStep * RShock_cvec + MShock_cvec.*uShock_cvec;
        RHSShock_cvec(1)    = 1;
        RHSShock_cvec(end)  = 0;


        MShock_cvec(1,1)     = 1;
        MShock_cvec(end,end) = 1;

        uSolutionShock_cvec  = RHSShock_cvec./MShock_cvec;
        uShock_cvec          = uSolutionShock_cvec;


        %% Expantion RHS and BC
        RHSExpantion_cvec       = -timeStep * RExpantion_cvec + MExpantion_cvec.*uExpantion_cvec;
        RHSExpantion_cvec(1)    = 0;
        RHSExpantion_cvec(end)  = 1;


        MExpantion_cvec(1,1)     = 1;
        MExpantion_cvec(end,end) = 1;

        uSolutionExpantion_cvec  = RHSExpantion_cvec./MExpantion_cvec;
        uExpantion_cvec          = uSolutionExpantion_cvec;

        %% Plotting Code
        if any(abs(Time - [0,0.25,0.5,0.75,1]) < 1e-6)
            legendOpts = {'NumColumns', 1, 'FontSize', 5, 'Box', 'off'};
            subplot(5,1,1)
            plot(x_vec',uSolution_cvec,'DisplayName', [num2str(Time) ' sec'])
            title(['Solution Using Full Mass Matrix Expresstion at $$\mu$$ = ' num2str(mu_vec(j))],'Interpreter','latex')
            xlabel('X','interpreter','latex')
            ylabel('U','interpreter','latex')
            hold on
            grid on
            axis equal
            legend('show', legendOpts{:});
            legend('Location', 'northwest');
    
            subplot(5,1,2)
            plot(x_vec',uSolutionLumped_mat,'DisplayName', [num2str(Time) ' sec'])
            title(['Solution Using Mass Lumping at $$\mu$$ = ' num2str(mu_vec(j))],'Interpreter','latex')
            % xlabel('X','interpreter','latex')
            ylabel('U','interpreter','latex')
            hold on
            grid on
            axis equal
            legend('show', legendOpts{:});
            legend('Location', 'northwest');
    
            subplot(5,1,3)
            plot(x_vec',uSolutionLumpedcvec_cvec,'DisplayName', [num2str(Time) ' sec'])
            title(['Solution Using Vector Mass Lumbing at $$\mu$$ = ' num2str(mu_vec(j))],'Interpreter','latex')
            % xlabel('X','interpreter','latex')
            ylabel('U','interpreter','latex')
            hold on
            grid on
            axis equal
            legend('show', legendOpts{:});
            legend('Location', 'northwest');
    
            subplot(5,1,4)
            plot(x_vec',uSolutionShock_cvec,'DisplayName', [num2str(Time) ' sec'])
            title(['Riemann Shock problem Solution at $$\mu$$ = ' num2str(mu_vec(j))],'Interpreter','latex')
            % xlabel('X','interpreter','latex')
            ylabel('U','interpreter','latex')
            hold on
            grid on
            axis equal
            legend('show', legendOpts{:});
            legend('Location', 'northwest');
    
            subplot(5,1,5)
            plot(x_vec',uSolutionExpantion_cvec,'DisplayName', [num2str(Time) ' sec'])
            title(['Riemann Expantion problem Solution at $$\mu$$ = ' num2str(mu_vec(j))],'Interpreter','latex')
            xlabel('X','interpreter','latex')
            ylabel('U','interpreter','latex')
            hold on
            grid on
            axis equal
            legend('show', legendOpts{:});
            legend('Location', 'northwest');
        end
        disp(Time)
        Time = Time + timeStep;
        itr = itr+1;
        disp(itr)
    end
end