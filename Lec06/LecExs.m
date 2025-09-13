clc;clearvars;close all;

%% Example 1

r1 = linspace(-0.5,0.5,11);
t1 = linspace(0,2,20);

figure;
for i = 1:length(r1)
    x = tanh(r1(i)).*t1 + r1(i);
    plot(x,t1,'DisplayName',['r = ' num2str(r1(i))])
    hold on
end
grid on
legend show
xlabel('x-axis','Interpreter','latex')
ylabel('t-axis','Interpreter','latex')
title('Curves for different r values','Interpreter','latex')

r = linspace(-10,10,100);
t = linspace(0,10,20);

figure;
for i = 1:length(t)

    w = tanh(r);
    x = tanh(r).*t(i) + r;
    plot3(x, t(i)*ones(size(r)), w);
    hold on
    grid on
    set(gca,'YDir','reverse')

end

xlabel('x-axis','Interpreter','latex')
ylabel('t-axis','Interpreter','latex')
zlabel('w-axis','Interpreter','latex')
title('3D Stepping Plot of $$w $$','Interpreter','latex')

% surface plotting
[R,T] = meshgrid(r,t);

W = tanh(R);
X = tanh(R).*T + R;

figure;
surf(X, T, W)

shading interp
set(gca,'YDir','reverse')

xlabel('x-axis','Interpreter','latex')
ylabel('t-axis','Interpreter','latex')
zlabel('w-axis','Interpreter','latex')
title('3D Surface Plot of w','Interpreter','latex')


%% Example 2
r1 = linspace(-2.5,0.5,11);
t1 = linspace(0,2,20);

figure;
for i = 1:length(r1)
    x = t1./(1+r1(i)^2) + r1(i);
    plot(x,t1,'DisplayName',['r = ' num2str(r1(i))])
    hold on
end
grid on
legend show
xlabel('x-axis','Interpreter','latex')
ylabel('t-axis','Interpreter','latex')
title('Curves for different r values','Interpreter','latex')

r = linspace(-10,10,100);
t = linspace(0,10,20);

figure;
for i = 1:length(t)

    w = 1./(1+r.^2);
    x = t(i)./(1+r.^2) + r;
    plot3(x, t(i)*ones(size(r)), w);
    hold on
    grid on
    set(gca,'YDir','reverse')

end

xlabel('x-axis','Interpreter','latex')
ylabel('t-axis','Interpreter','latex')
zlabel('w-axis','Interpreter','latex')
title('3D Stepping Plot of $$w $$','Interpreter','latex')

% surface plotting
[R,T] = meshgrid(r,t);

W = 1./(1+R.^2);
X = T./(1+R.^2) + R;

figure;
surf(X, T, W)

shading interp
set(gca,'YDir','reverse')

xlabel('x-axis','Interpreter','latex')
ylabel('t-axis','Interpreter','latex')
zlabel('w-axis','Interpreter','latex')
title('3D Surface Plot of w','Interpreter','latex')


