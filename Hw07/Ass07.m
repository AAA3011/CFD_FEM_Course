clc;clearvars;close all;

domainLength    = 1;
Pe_vec          = [1 10 100 500];
numElements_vec = [10 100];
numRows         = length(Pe_vec);

for j = 1:length(numElements_vec)
    figure;
    for i = 1:length(Pe_vec)

    elementLength = domainLength/numElements_vec(j);
    x_vec         = 0:elementLength:1;
    %% Solution Using Leat-Square Methode
    FLS   = 0;
    ThoLS = 1;
    TLSSolution_cvec = methode(numElements_vec(j),elementLength,Pe_vec(i),FLS,ThoLS);

    %% Solution Using Standard Galarkin Methode
    FSG   = 1;
    ThoSG = 0;
    TSGSolution_cvec = methode(numElements_vec(j),elementLength,Pe_vec(i),FSG,ThoSG);

    %% Solution Using Galarkin/Least-square Methode9
    FGLS   = 1;
    ThoGLS = 1/sqrt((2*Pe_vec(i)/elementLength)^2+(4*-1/elementLength)^2); 
    TGLSSolution_cvec = methode(numElements_vec(j),elementLength,Pe_vec(i),FGLS,ThoGLS);

    %% Exact Solution
    TExactSolution_vec = (exp(Pe_vec(i)*x_vec))/(exp(Pe_vec(i))-1);    
    %% Plotting
    subplot(numRows,1,i)
    plot(x_vec,TLSSolution_cvec,x_vec,TSGSolution_cvec,x_vec,TGLSSolution_cvec,x_vec,TExactSolution_vec)
    % plot(x_vec,TSGSolution_cvec,x_vec,TGLSSolution_cvec,x_vec,TExactSolution_vec)
    title('1D Steady Heat Conduction in a Rod with Upwinding','Interpreter','latex')
    title(['1D Steady Heat Conduction in a Rod with Upwinding at $Pe = ',num2str(Pe_vec(i)),'$',' and $ Number of Elements = ',num2str(numElements_vec(j)),'$'], 'Interpreter', 'latex');
    xlabel('X ($$m$$)','Interpreter','latex')
    ylabel('Temperature ($$K$$)','Interpreter','latex')
    grid on
    legend('LS','SG','GLS','Exact','Location','northwest')
    % legend('SG','GLS','Exact','Location','northwest')

    end
end

%% Clustered Case
numElements  = 10;
rate         = 0.94;
domainLength = 1;
PeC_vec      = [1 10 100 500];
numRowsC     = length(PeC_vec);
plotNum      = 1;

figure;
for i = 1:length(PeC_vec)

    %% Solution Using Standard Galarkin Methode
    FSG   = 1;
    ThoSG = 0;
    [TSGSolutionC_cvec,PeLocal_vec,x_vec] = ClusteredCase(numElements,domainLength,rate,FSG,ThoSG,PeC_vec(i));
    
    %% Solution Using Galarkin/Least-square Methode9
    FGLS   = 1;
    ThoGLS = 2;  % This is not the actual values however in the ClusteredCase function, the if condition ensures the correct value is calculated
    [TGLSSolutionC_cvec,~,~] = ClusteredCase(numElements,domainLength,rate,FGLS,ThoGLS,PeC_vec(i));
    
    %% Exact Solution
    TExactSolutionC_vec = (exp(PeC_vec(i)*x_vec))/(exp(PeC_vec(i))-1);

    %% Plotting
    subplot(numRowsC,2,plotNum)
    plot(x_vec,TSGSolutionC_cvec,x_vec,TGLSSolutionC_cvec,x_vec,TExactSolutionC_vec)
    title(['1D Steady Heat Conduction in a Rod with Upwinding at $Pe = ',num2str(PeC_vec(i)),'$',' and $ Number of Elements = ',num2str(numElements),'$'], 'Interpreter', 'latex');
    xlabel('X ($$m$$)','Interpreter','latex')
    ylabel('Temperature ($$K$$)','Interpreter','latex')
    grid on
    legend('SG','GLS','Exact','Location','northwest')

    subplot(numRowsC,2,plotNum+1)
    scatter(x_vec,PeLocal_vec)
    hold on
    plot(x_vec,PeLocal_vec)
    hold off
    grid on
    title(['local Peclet number at $Pe = ',num2str(PeC_vec(i)),'$',' and $ Number of Elements = ',num2str(numElements),'$'], 'Interpreter', 'latex');
    xlabel('X ($$m$$)','Interpreter','latex')
    ylabel('$$Pe_{local}$$','Interpreter','latex')
    plotNum = plotNum+2;
end
