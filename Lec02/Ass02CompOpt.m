clc;close all;clearvars;
tic

%% Connectivity Matrix 

numElements = 100; 
numNodes    = numElements + 1;
connectivityMatrix_mat = zeros(numElements,2);

connectivityMatrix_mat(:,1) = 1:numElements ;
connectivityMatrix_mat(:,2) = connectivityMatrix_mat(:,1)+1;

node1_cvec = connectivityMatrix_mat(:,1);
node2_cvec = connectivityMatrix_mat(:,2);
    
%% Globla Matrix
L = 1;
totalLength   = 10;

indexAreaHalf = floor(numNodes/2);
indexKHalf    = floor(numElements/2);

globlaMatrix_mat = zeros(numNodes,numNodes);
area_vec         = zeros(1,numNodes);
K_vec            = zeros(1,numElements);

elementLength   = totalLength / numElements; 
X_vec           = 0:elementLength:totalLength;
rho_vec         = ones(1,numNodes)*1.225;  
phiXL_vec       = [500 1500 2500 3000 3400 3430 3500 500];   % the last value was put to be 500 so that the last curves are the same that in the assienment
mach_mat        = zeros(numNodes,length(phiXL_vec));
temperature_mat = zeros(numNodes,length(phiXL_vec));
pressure_mat    = zeros(numNodes,length(phiXL_vec));

for i = 1 : length(phiXL_vec)
    while (L > 10^-6)
        
        K_vec(1:indexKHalf) = ((0.1*rho_vec(1:indexKHalf)*totalLength/elementLength) - (0.15 * rho_vec(1:indexKHalf) * totalLength^2 / (6*elementLength^2)) .* ( (1- 2*(X_vec(1:indexKHalf)+elementLength)/totalLength).^3  - (1- 2*X_vec(1:indexKHalf)/totalLength).^3));
        area_vec(1:indexAreaHalf) = 0.1 * totalLength + 0.15 * totalLength * (1- 2 * X_vec(1:indexAreaHalf)/totalLength).^2 ;
        
        K_vec(indexKHalf+1:numElements) = ((0.1*rho_vec(indexKHalf+1:numElements)*totalLength/elementLength) + (0.05 * rho_vec(indexKHalf+1:numElements) * totalLength^2 / (6*elementLength^2)) .* ( (2*(X_vec(indexKHalf+1:numElements)+elementLength)/totalLength - 1 ).^3  - (2*X_vec(indexKHalf+1:numElements)/totalLength  - 1).^3));
        area_vec(indexAreaHalf+1:numNodes) = 0.1 * totalLength + 0.05 * totalLength * (2 * X_vec(indexAreaHalf+1:numNodes)/totalLength - 1).^2 ;
    
    
        rows_cvec = [node1_cvec;node2_cvec;node1_cvec;node2_cvec];
        cols_cvec = [node1_cvec;node2_cvec;node2_cvec;node1_cvec];
        vals_cvec = [K_vec(:);K_vec(:);-K_vec(:);-K_vec(:)];
        
        globlaMatrix_mat = sparse(rows_cvec, cols_cvec, vals_cvec, numNodes, numNodes);
        globlaMatrix_mat = full(globlaMatrix_mat);
    
    
        %% Boundary Conditions
        
        phiX0 = 0 ;
    
        RHS_cvec  = zeros(numNodes,1);
        
        RHS_cvec(end) = RHS_cvec(end) +  phiXL_vec(i);
        RHS_cvec(1)   = RHS_cvec(1)   +  phiX0;
        
        globlaMatrix_mat(1,:)     = 0;
        globlaMatrix_mat(1,1)     = 1;
        globlaMatrix_mat(end,:)   = 0;
        globlaMatrix_mat(end,end) = 1;
    
        phiConstants_cvec = globlaMatrix_mat \ RHS_cvec ; 
            
        % Velocity Using FEM
    
        uKVal_vec  =  ones(1,numElements) * (elementLength/3) ;
        uvals_cvec = [uKVal_vec(:);uKVal_vec(:);uKVal_vec(:)/2;uKVal_vec(:)/2];
        
        uGloblaMatrix_mat = sparse(rows_cvec, cols_cvec, uvals_cvec, numNodes, numNodes);
        uGloblaMatrix_mat = full(uGloblaMatrix_mat);
        
        uRHS_cvec = zeros(numNodes,1) ;
        
        uRHS_cvec(1:numElements) = uRHS_cvec(1:numElements) + (phiConstants_cvec(2:numNodes) - phiConstants_cvec(1:numElements)) / 2 ;
        uRHS_cvec(2:numNodes)    = uRHS_cvec(2:numNodes)    + (phiConstants_cvec(2:numNodes) - phiConstants_cvec(1:numElements)) / 2 ;
        
        velocityFEM_cvec = uGloblaMatrix_mat \ uRHS_cvec;
        
        %% rho Calculations
    
        gasConstant = 287;
        pressure0   = 3*10^5 ; 
        gamma       = 1.4;
        rho0        = 1.225;
    
        temperature0  = pressure0/(gasConstant * rho0);
        airspeed0     = sqrt(gamma * gasConstant * temperature0) ;
        mach0_vec     = velocityFEM_cvec' ./ airspeed0;
        rhoNew_vec    = rho0 * (1-(((gamma-1)./2).*mach0_vec.^2)) .^ (1/(gamma-1));
        
        %% Condition
        
        R_vec   = 2*(rhoNew_vec - rho_vec)./(rhoNew_vec + rho_vec);
        L       = sum(R_vec.^2);
        rho_vec = rhoNew_vec ;
    
    end
    L = 1;
    mach_mat(:,i)        = sqrt((-2*mach0_vec.^2)./((mach0_vec.^2 * (gamma -1))-2));
    temperature_mat(:,i) = temperature0 * (1+((gamma-1)/2).*(mach_mat(:,i).^2)).^-1 ;
    pressure_mat(:,i)    = rho_vec' .* gasConstant .* temperature_mat(:,i) ;
    
end

%% Plotting Code

figure(1)
plot(X_vec,phiConstants_cvec)
title('Potential Function','interpreter','latex')
ylabel('$$\phi(x)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on

figure(2)
subplot(4,1,1)
plot(X_vec,area_vec,'k','LineWidth', 4)
title('Convergent Divergent Nozzel shape','interpreter','latex')
ylabel('Area $$(m^2)$$','interpreter','latex')
ylim([0 max(area_vec)* 1.3])
grid on

subplot(4,1,2)
plot(X_vec,rho_vec)
title('Density Variation along the Nozzle','interpreter','latex')
ylabel('$$\rho \ (Kg/m^3)$$','interpreter','latex')
grid on


subplot(4,1,3)
plot(X_vec,temperature_mat(:,end))
title('Temperature Variation along the Nozzle','interpreter','latex')
ylabel('T $$(K)$$','interpreter','latex')
grid on

subplot(4,1,4)
plot(X_vec,pressure_mat(:,end)./10^5)
title('Pressure Variation along the Nozzle','interpreter','latex')
ylabel('P $$(bar)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on

figure(3)
plot(X_vec,mach_mat(:,1:end-1))
title('Mach number Variation for Different Potential Functions at exit','interpreter','latex')
ylabel('$$M$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
legendStrings = arrayfun(@(phi) sprintf('$$\\phi_{exit}(x) = %d  $$',phi),phiXL_vec(1,1:end-1), 'UniformOutput', false);
legend(legendStrings, 'Interpreter','latex')
grid on

figure(4)
plot(X_vec,velocityFEM_cvec)
title('Velocity Variation Using (SG FEM)','interpreter','latex')
ylabel('Velocity $$(m/s)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on
toc