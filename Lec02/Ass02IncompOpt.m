tic
clc;close all;clearvars;

%% Connectivity Matrix 

numElements = 100 ; 
numNodes    = numElements + 1;
connectivityMatrix_mat = zeros(numElements,2);

connectivityMatrix_mat(:,1) = 1:numElements ;
connectivityMatrix_mat(:,2) = connectivityMatrix_mat(:,1)+1;

%% Area and K calculations

rho = 1.225;  
totalLength   = 10;
elementLength = totalLength / numElements; 
 
X_vec  = 0:elementLength:totalLength;

indexAreaHalf = floor(numNodes/2);
indexKHalf    = floor(numElements/2);
area_vec      = zeros(1,numNodes);
K_vec         = zeros(1,numElements);

K_vec(1:indexKHalf) = ((0.1*rho*totalLength/elementLength) - (0.15 * rho * totalLength^2 / (6*elementLength^2)) * ( (1- 2*(X_vec(1:indexKHalf)+elementLength)/totalLength).^3  - (1- 2*X_vec(1:indexKHalf)/totalLength).^3));
area_vec(1:indexAreaHalf) = 0.1 * totalLength + 0.15 * totalLength * (1- 2 * X_vec(1:indexAreaHalf)/totalLength).^2 ;

K_vec(indexKHalf+1:numElements) = ((0.1*rho*totalLength/elementLength) + (0.05 * rho * totalLength^2 / (6*elementLength^2)) * ( (2*(X_vec(indexKHalf+1:numElements)+elementLength)/totalLength - 1 ).^3  - (2*X_vec(indexKHalf+1:numElements)/totalLength  - 1).^3));
area_vec(indexAreaHalf+1:numNodes) = 0.1 * totalLength + 0.05 * totalLength * (2 * X_vec(indexAreaHalf+1:numNodes)/totalLength - 1).^2 ;

%% Global Matrix

node1_cvec = connectivityMatrix_mat(:,1);
node2_cvec = connectivityMatrix_mat(:,2);

rows_cvec = [node1_cvec;node2_cvec;node1_cvec;node2_cvec];
cols_cvec = [node1_cvec;node2_cvec;node2_cvec;node1_cvec];
vals_cvec = [K_vec(:);K_vec(:);-K_vec(:);-K_vec(:)];

globlaMatrix_mat = sparse(rows_cvec, cols_cvec, vals_cvec, numNodes, numNodes);
globlaMatrix_mat = full(globlaMatrix_mat);

%% Boundary Conditions

pressureXL_vec = [2.95 2.9 2.85 2.8 2.7 2.97]* 10^5;
pressure0      = 3 * 10^5;
phiX0          = 0;

phiDashXL_vec = sqrt(2* (pressure0 - pressureXL_vec) ./ rho);

RHS_mat  = zeros(numNodes,length(pressureXL_vec));

RHS_mat(end,:) = RHS_mat(end,:) + rho * area_vec(end) .* phiDashXL_vec;
RHS_mat(1,:)   = RHS_mat(1,:)   +  phiX0;

globlaMatrix_mat(1,:) = 0;
globlaMatrix_mat(1,1) = 1;
    
phiConstants_mat = globlaMatrix_mat \ RHS_mat(:,:) ; 
    

%% Requierments 

% Velocity Using FDM
velocityFDM_mat                  = zeros(numNodes,length(pressureXL_vec));
velocityFDM_mat(1:numElements,:) = (phiConstants_mat(2:numNodes,:) - phiConstants_mat(1:numElements,:)) ./ elementLength;
velocityFDM_mat(end,:)           = phiDashXL_vec;

% Velocity Using FEM
uKVal_vec  =  ones(1,numElements) * (elementLength/3) ;
uvals_cvec = [uKVal_vec(:);uKVal_vec(:);uKVal_vec(:)/2;uKVal_vec(:)/2];

uGloblaMatrix_mat = sparse(rows_cvec, cols_cvec, uvals_cvec, numNodes, numNodes);
uGloblaMatrix_mat = full(uGloblaMatrix_mat);

uRHS_mat = zeros(numNodes,length(pressureXL_vec)) ;

uRHS_mat(1:numElements,:) = uRHS_mat(1:numElements,:) + (phiConstants_mat(2:numNodes,:) - phiConstants_mat(1:numElements,:)) / 2 ;
uRHS_mat(2:numNodes,:)    = uRHS_mat(2:numNodes,:)    + (phiConstants_mat(2:numNodes,:) - phiConstants_mat(1:numElements,:)) / 2 ;

velocityFEM_mat = uGloblaMatrix_mat \ uRHS_mat;

% Velocity Exact Solution
massFlowRateXL_vec = phiDashXL_vec * rho * area_vec(end);  
velocityExact_mat  = massFlowRateXL_vec ./ (rho * area_vec');

% Mass Flow Rate Calculation
massFlowRateFDM_mat   = rho .* velocityFDM_mat   .* area_vec' ;
massFlowRateFEM_mat   = rho .* velocityFEM_mat   .* area_vec' ;
massFlowRateExact_mat = rho .* velocityExact_mat .* area_vec' ;

% Flow Parameters
gamma          = 1.4;
gasConst       = 287 ;

pressure_mat     = pressure0 - 0.5 * rho * velocityFEM_mat.^2 ;
temperature_mat  = pressure_mat  ./ (rho * gasConst) ;
airSpeed_mat     = sqrt(gamma * gasConst * temperature_mat);
mach_mat         = velocityFEM_mat ./ airSpeed_mat;

%% Plotting Code

figure(1)
plot(X_vec,phiConstants_mat(:,end))
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
plot(X_vec,temperature_mat(:,end))
title('Temperature Variation along the Nozzle','interpreter','latex')
ylabel('T $$(K)$$','interpreter','latex')
grid on

subplot(3,1,3)
plot(X_vec,pressure_mat(:,end)/10^5)
title('Pressure Variation along the Nozzle','interpreter','latex')
ylabel('P $$(bar)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on

figure(3)
subplot(1,3,1)
plot(X_vec,velocityFDM_mat(:,end))
title('Velocity Variation Using (Foward FDM)','interpreter','latex')
ylabel('Velocity $$(m/s)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on

subplot(1,3,2)
plot(X_vec,velocityFEM_mat(:,end))
title('Velocity Variation Using (SG FEM)','interpreter','latex')
ylabel('Velocity $$(m/s)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on

subplot(1,3,3)
plot(X_vec,velocityFEM_mat(:,end),X_vec,velocityFDM_mat(:,end),X_vec,velocityExact_mat(:,end),'o','MarkerSize',4)
title('Comparison Between Velocity variation Using Different Techniques','interpreter','latex')
ylabel('Velocity $$(m/s)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
legend('FEM','FDM','Exact','interpreter','latex')
grid on

figure(4)
plot(X_vec,massFlowRateFDM_mat(:,end),X_vec,massFlowRateFEM_mat(:,end),X_vec,massFlowRateExact_mat(:,end))
title('Comparison Between Mass Flow Rate Using Different Techniques','interpreter','latex')
ylabel('$$\dot{m}$$ $$(kg/s)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
legend('FDM','FEM','Exact','interpreter','latex')
ylim([min(massFlowRateFDM_mat(:,end))*0.7 max(massFlowRateFDM_mat(:,end))*1.3 ])
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