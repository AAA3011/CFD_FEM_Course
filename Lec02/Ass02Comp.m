clc;close all;clearvars;
tic

%% Connectivity Matrix 

numElements = 100; 
numNodes    = numElements + 1;
connectivityMatrix_mat = zeros(numElements,2);

for i = 1 : numElements

    connectivityMatrix_mat(i,1) = i;
    connectivityMatrix_mat(i,2) = i+1;

end

%% Globla Matrix
rho_vec = ones(1,numNodes)*1.225;  
totalLength   = 10;
elementLength = totalLength / numElements ; 
X_vec         = zeros(1,numNodes);

for i = 1:numElements

    X_vec(i+1) = X_vec(i) + elementLength;

end

area_vec = zeros(1,numNodes);
K_vec    = zeros(1,numElements);

L = 1;

phiXL_vec       = [500 1500 2500 3000 3400 3430 3500 500];  % the last value was put to be 500 so that the last curves are the same that in the assienment
mach_mat        = zeros(numNodes,length(phiXL_vec));
temperature_mat = zeros(numNodes,length(phiXL_vec));
pressure_mat    = zeros(numNodes,length(phiXL_vec));

for j = 1:length(phiXL_vec)
    while (L > 10^-6)
        globlaMatrix_mat  = zeros(numNodes,numNodes);
        RHS_cvec          = zeros(numNodes,1);
        uGlobalMatrix_mat = zeros(numNodes,numNodes);
        uRHS_cvec         = zeros(numNodes,1);

        for i = 1:numElements 

            if (0 <= X_vec(i) && X_vec(i) <= (totalLength /2))
                
                K_vec(i)    = ((0.1*rho_vec(i)*totalLength/elementLength) - (0.15 * rho_vec(i) * totalLength^2 / (6*elementLength^2)) * ( (1- 2*X_vec(i+1)/totalLength)^3  - (1- 2*X_vec(i)/totalLength)^3));
                area_vec(i) = 0.1 * totalLength + 0.15 * totalLength * (1- 2 * X_vec(i)/totalLength)^2 ;
        
            elseif ((totalLength /2) < X_vec(i) && X_vec(i) <= totalLength) 
        
                K_vec(i)    = ((0.1*rho_vec(i)*totalLength/elementLength) + (0.05 * rho_vec(i) * totalLength^2 / (6*elementLength^2)) * ( (2*X_vec(i+1)/totalLength - 1 )^3  - (2*X_vec(i)/totalLength  - 1)^3));
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
        
        RHS_cvec(end) =  phiXL_vec(j);
        RHS_cvec(1)   =  phiX0;
            
        globlaMatrix_mat(1,:)     = 0;
        globlaMatrix_mat(end,:)   = 0;
        globlaMatrix_mat(1,1)     = 1;
        globlaMatrix_mat(end,end) = 1;
        
        phiConstants_cvec = globlaMatrix_mat \ RHS_cvec ; 
            
        % Velocity Using FEM  

        for i = 1:numElements
            
            uGlobalMatrix_mat(i+1,i+1) = uGlobalMatrix_mat(i+1,i+1) + elementLength / 3;
            uGlobalMatrix_mat(i,i)     = uGlobalMatrix_mat(i,i)     + elementLength / 3;
            uGlobalMatrix_mat(i+1,i)   = uGlobalMatrix_mat(i+1,i)   + elementLength / 6;
            uGlobalMatrix_mat(i,i+1)   = uGlobalMatrix_mat(i,i+1)   + elementLength / 6;
                
            uRHS_cvec(i)   = uRHS_cvec(i)   + (phiConstants_cvec(i+1) - phiConstants_cvec(i)) / 2 ;
            uRHS_cvec(i+1) = uRHS_cvec(i+1) + (phiConstants_cvec(i+1) - phiConstants_cvec(i)) / 2 ;
                
        end
            
        velocityFEM_cvec = uGlobalMatrix_mat \ uRHS_cvec;
        
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
    %% Flow parameters
    L = 1;
    mach_mat(:,j)        = sqrt((-2*mach0_vec.^2)./((mach0_vec.^2 * (gamma -1))-2));
    temperature_mat(:,j) = temperature0 * (1+((gamma-1)/2).*(mach_mat(:,j).^2)).^-1 ;
    pressure_mat(:,j)    = rho_vec' .* gasConstant .* temperature_mat(:,j);

end

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
legendStrings = arrayfun(@(phi) sprintf('$$\\phi(x) = %d  $$',phi),phiXL_vec(1,1:end-1), 'UniformOutput', false);
legend(legendStrings, 'Interpreter','latex')
grid on

figure(4)
plot(X_vec,velocityFEM_cvec)
title('Velocity Variation Using (SG FEM)','interpreter','latex')
ylabel('Velocity $$(m/s)$$','interpreter','latex')
xlabel('X $$(m)$$','interpreter','latex')
grid on

toc