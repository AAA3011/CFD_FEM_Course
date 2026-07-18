% usefulPlots - Useful Plots.
%
% FILE: usefulPlots.m
% DESCRIPTION:
% Collection of small plotting helper scripts that load precomputed norm
% histories and produce diagnostic figures for convergence and iteration
% counts across time steps.
%
% INPUTS:
% Expects MAT files like `norms1.mat`, `norms500.mat`, `stepNorms.mat` in cwd
%
% OUTPUTS:
% On-screen figures and saved images (if saving enabled inside script)
% % First Time Step iterations Convergence
%
% Syntax: usefulPlots()
%
% Inputs: none
%
% Outputs: none
%

data = load('norms1.mat');
data = data.normsOut;

figure;
semilogy(1:size(data,1),data(:,1),'-o', 'LineWidth', 1.2, 'MarkerSize', 4);hold on;
semilogy(1:size(data,1),data(:,2),'-s', 'LineWidth', 1.2, 'MarkerSize', 4)
semilogy(1:size(data,1),data(:,3),'-^', 'LineWidth', 1.2, 'MarkerSize', 4);hold off;

grid on;
xlabel('Time step', 'Interpreter', 'latex');
ylabel('Converged norm value', 'Interpreter', 'latex');
title(['Convergence History at time step = ',num2str(1) ], 'Interpreter', 'latex');
legend({'$||du||_2$', '$||dv||_2$', '$||dp||_2$'},'Interpreter', 'latex', 'Location', 'best');

%% 500 Time Step iterations Convergence
data = load('norms500.mat');
data = data.normsOut;

figure;
semilogy(1:size(data,1),data(:,1),'-o', 'LineWidth', 1.2, 'MarkerSize', 4);hold on;
semilogy(1:size(data,1),data(:,2),'-s', 'LineWidth', 1.2, 'MarkerSize', 4)
semilogy(1:size(data,1),data(:,3),'-^', 'LineWidth', 1.2, 'MarkerSize', 4);hold off;

grid on;
xlabel('Time step', 'Interpreter', 'latex');
ylabel('Converged norm value', 'Interpreter', 'latex');
title(['Convergence History at time step = ',num2str(500) ], 'Interpreter', 'latex');
legend({'$||du||_2$', '$||dv||_2$', '$||dp||_2$'},'Interpreter', 'latex', 'Location', 'best');

%% 1000 Time Step iterations Convergence
data = load('norms1000.mat');
data = data.normsOut;

figure;
semilogy(1:size(data,1),data(:,1),'-o', 'LineWidth', 1.2, 'MarkerSize', 4);hold on;
semilogy(1:size(data,1),data(:,2),'-s', 'LineWidth', 1.2, 'MarkerSize', 4)
semilogy(1:size(data,1),data(:,3),'-^', 'LineWidth', 1.2, 'MarkerSize', 4);hold off;

grid on;
xlabel('Time step', 'Interpreter', 'latex');
ylabel('Converged norm value', 'Interpreter', 'latex');
title(['Convergence History at time step = ',num2str(1000) ], 'Interpreter', 'latex');
legend({'$||du||_2$', '$||dv||_2$', '$||dp||_2$'},'Interpreter', 'latex', 'Location', 'best');


%% stepNorms Plotting
data = load('stepNorms.mat');
data = data.stepNorms_mat;

figure;
semilogy(1:size(data,1),data(:,1),'-o', 'LineWidth', 1.2, 'MarkerSize', 4);hold on;
semilogy(1:size(data,1),data(:,2),'-s', 'LineWidth', 1.2, 'MarkerSize', 4)
semilogy(1:size(data,1),data(:,3),'-^', 'LineWidth', 1.2, 'MarkerSize', 4);hold off;

grid on;
xlabel('Time step', 'Interpreter', 'latex');
ylabel('Converged norm value', 'Interpreter', 'latex');
title('Convergence History', 'Interpreter', 'latex');
legend({'$||du||_2$', '$||dv||_2$', '$||dp||_2$'},'Interpreter', 'latex', 'Location', 'best');

%% line Convergence History
figure;
semilogy(1:size(data,1),data(:,4), 'r-o', 'LineWidth', 1.2, 'MarkerSize', 3)
xlabel('Time','Interpreter','latex');
ylabel('$||RHS||_2$ at convergence','Interpreter','latex');
title('Residual per Time Step','Interpreter','latex');
grid on;

%% Number of Iterations per time Steps
figure;
plot(1:size(data,1),data(:,5), 'k-o', 'LineWidth', 1.2, 'MarkerSize', 3);
xlabel('Time','Interpreter','latex');
ylabel('Newton iterations','Interpreter','latex');
title('Newton Iterations per Time Step','Interpreter','latex');
grid on;

%% Divergence History per time Steps
figure;
plot(1:size(data,1),data(:,6), '-o', 'LineWidth', 0.7, 'MarkerSize', 3);
xlabel('Time','Interpreter','latex');
ylabel('Time Steps','Interpreter','latex');
title('Divergence History per Time Step','Interpreter','latex');
grid on;
