tic
clc;close all;clearvars;

%% Connectivity Matrix 

numElements = 100 ; 
numNodes    = numElements + 1;
connectivityMatrix_mat = zeros(numElements,2);

for i = 1 : numElements

    connectivityMatrix_mat(i,1) = i;
    connectivityMatrix_mat(i,2) = i+1;

end

%% Globla Matrix

rho = 1.225;  
totalLength   = 10;
elementLength = totalLength / numElements ; 
X_vec         = zeros(1,numNodes);

for i = 1:numElements

    X_vec(i+1) = X_vec(i) + elementLength;

end

area_vec = zeros(1,numNodes);
K_vec    = zeros(1,numElements);
globlaMatrix_mat = zeros(numNodes,numNodes);

for i = 1:numElements

    if (0 <= X_vec(i) && X_vec(i) <= (totalLength /2))
        
        K_vec(i) = ((0.1*rho*totalLength/elementLength) - (0.15 * rho * totalLength^2 / (6*elementLength^2)) * ( (1- 2*X_vec(i+1)/totalLength)^3  - (1- 2*X_vec(i)/totalLength)^3));
        area_vec(i) = 0.1 * totalLength + 0.15 * totalLength * (1- 2 * X_vec(i)/totalLength)^2 ;

    elseif ((totalLength /2) < X_vec(i) && X_vec(i) <= totalLength) 

        K_vec(i) = ((0.1*rho*totalLength/elementLength) + (0.05 * rho * totalLength^2 / (6*elementLength^2)) * ( (2*X_vec(i+1)/totalLength - 1 )^3  - (2*X_vec(i)/totalLength  - 1)^3));
        area_vec(i) = 0.1 * totalLength + 0.05 * totalLength * (2 * X_vec(i)/totalLength - 1)^2 ;

    end
        globlaMatrix_mat(i,i)     = globlaMatrix_mat(i,i)     + K_vec(i);
        globlaMatrix_mat(i+1,i+1) = globlaMatrix_mat(i+1,i+1) + K_vec(i);
        globlaMatrix_mat(i,i+1)   = globlaMatrix_mat(i,i+1)   - K_vec(i);
        globlaMatrix_mat(i+1,i)   = globlaMatrix_mat(i+1,i)   - K_vec(i);

end
area_vec(end) = 0.1 * totalLength + 0.05 * totalLength * (2 * X_vec(end)/totalLength - 1)^2 ;

%% Boundary Conditions

phiX0 = 0 ;
pressureXL_vec = [2.95 2.9 2.85 2.8 2.7 2.97 ]* 10^5 ;
pressure0  = 3 * 10^5;

phiDashXL_vec = sqrt(2* (pressure0 - pressureXL_vec) ./ rho);
mach_mat      = zeros(numNodes,length(pressureXL_vec));
for j = 1:length(pressureXL_vec)

    RHS_cvec  = zeros(numNodes,1);
    RHS_cvec(end) = RHS_cvec(end) + rho * area_vec(end) * phiDashXL_vec(j);
    RHS_cvec(1)   = RHS_cvec(1)   + phiX0;
    
    globlaMatrix_mat(1,:) = 0;
    globlaMatrix_mat(1,1) = 1;
    
    phiConstants_cvec = globlaMatrix_mat \ RHS_cvec ; 
    
    
    %% Requierments 
    
    % Velocity Using FDM
    velocityFDM_cvec = zeros(numNodes,1);
    for i = 1:numElements
    
        velocityFDM_cvec(i) = (phiConstants_cvec(i+1) - phiConstants_cvec(i)) / elementLength;
    
    end
    velocityFDM_cvec(end) = phiDashXL_vec(j);
    
    % Velocity Using FEM
    uGlobalMatrix_mat = zeros(numNodes,numNodes);
    uRHS_cvec         = zeros(numNodes,1);
    for i = 1:numElements
    
        uGlobalMatrix_mat(i+1,i+1) = uGlobalMatrix_mat(i+1,i+1) + elementLength / 3;
        uGlobalMatrix_mat(i,i)     = uGlobalMatrix_mat(i,i)     + elementLength / 3;
        uGlobalMatrix_mat(i+1,i)   = uGlobalMatrix_mat(i+1,i)   + elementLength / 6;
        uGlobalMatrix_mat(i,i+1)   = uGlobalMatrix_mat(i,i+1)   + elementLength / 6;
    
        uRHS_cvec(i)   = uRHS_cvec(i) + (phiConstants_cvec(i+1) - phiConstants_cvec(i)) / 2 ;
        uRHS_cvec(i+1) = uRHS_cvec(i+1) + (phiConstants_cvec(i+1) - phiConstants_cvec(i)) / 2 ;
    
    end
    
    velocityFEM_cvec = uGlobalMatrix_mat \ uRHS_cvec;
    
    % Velocity Exact Solution
    massFlowRateXL = phiDashXL_vec(j) * rho * area_vec(end) ;  
    velocityExact_vec = massFlowRateXL ./ (rho * area_vec) ;
    
    % Mass Flow Rate Calculation
    massFlowRateFDM_vec = rho .* velocityFDM_cvec' .* area_vec ;
    massFlowRateFEM_vec = rho .* velocityFEM_cvec' .* area_vec ;
    massFlowRateExact_vec = rho .* velocityExact_vec .* area_vec ;

    % Flow Parameters
    gamma          = 1.4;
    gasConst       = 287 ;

    pressure_cvec    = pressure0 - 0.5 * rho * velocityFEM_cvec.^2 ;
    temperature_cvec = pressure_cvec  ./ (rho * gasConst) ;
    airSpeed_cvec    = sqrt(gamma * gasConst * temperature_cvec);
    mach_mat(:,j)    = velocityFEM_cvec ./ airSpeed_cvec;

end

%% Plotting Code

figure(1)
plot(X_vec,phiConstants_cvec)
title('Potential Function','interpreter','latex')
ylabel('$$\phi(x)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on

figure(2)
subplot(3,1,1)
plot(X_vec,area_vec,'k','LineWidth', 4)
title('Convergent Divergent Nozzle shape','interpreter','latex')
ylabel('Area $$(m^2)$$','interpreter','latex')
ylim([0 max(area_vec)* 1.3])
grid on

subplot(3,1,2)
plot(X_vec,temperature_cvec)
title('Temperature Variation along the Nozzle','interpreter','latex')
ylabel('T $$(K)$$','interpreter','latex')
grid on

subplot(3,1,3)
plot(X_vec,pressure_cvec/10^5)
title('Pressure Variation along the Nozzle','interpreter','latex')
ylabel('P $$(bar)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on

figure(3)
subplot(1,3,1)
plot(X_vec,velocityFDM_cvec)
title('Velocity Variation Using (Foward FDM)','interpreter','latex')
ylabel('Velocity $$(m/s)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on

subplot(1,3,2)
plot(X_vec,velocityFEM_cvec)
title('Velocity Variation Using (SG FEM)','interpreter','latex')
ylabel('Velocity $$(m/s)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on

subplot(1,3,3)
plot(X_vec,velocityFEM_cvec,X_vec,velocityFDM_cvec,X_vec,velocityExact_vec,'o','MarkerSize',4)
title('Comparison Between Velocity variation Using Different Techniques','interpreter','latex')
ylabel('Velocity $$(m/s)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
legend('FEM','FDM','Exact','interpreter','latex')
grid on

figure(4)
plot(X_vec,massFlowRateFDM_vec,X_vec,massFlowRateFEM_vec,X_vec,massFlowRateExact_vec)
title('Comparison Between Mass Flow Rate Using Different Techniques','interpreter','latex')
ylabel('$$\dot{m}$$ $$(kg/s)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
legend('FDM','FEM','Exact','interpreter','latex')
ylim([min(massFlowRateFDM_vec)*0.7 max(massFlowRateFDM_vec)*1.3 ])
grid on

figure(5)
plot(X_vec,mach_mat)
title('Mach number Variation for Different Back pressure Values','interpreter','latex')
ylabel('$$M$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
legendStrings = arrayfun(@(p) sprintf('$$p_b = %.2f \\, bar$$', p/1e5),pressureXL_vec, 'UniformOutput', false);
legend(legendStrings, 'Interpreter','latex')
grid on

toc