%% ELEC 518 Problem 6
%
% Ender Erkaya
%
% January 2022
%
clearvars
%% Settings
h    = 1e-2;
x0   = 1;

%% Function Handles
xorig = @(u) cos(u);
fun   = @(x,u) -100*(x-cos(u))-sin(u);
Jfun  = @(x) -100;

%% 6b
Tmax  = 1;
figure
[xAB2,tAB2] = AB2(fun, Jfun, x0, 1e-2, Tmax);
hold on
plot(tAB2,xAB2,'m','LineWidth',2);
hold on
plot(tAB2,xorig(tAB2),'b','LineWidth',2);
[xAB2,tAB2] = AB2(fun, Jfun, x0, 1.07e-2, Tmax);
hold on
plot(tAB2,xAB2,'r','LineWidth',2);
grid on
legend('h=1e-2','exact','h=1.07e-2')
xlabel('t','Interpreter','latex');
ylabel('x(t)','Interpreter','latex');
title('AB2 vs Exact Solution','Interpreter','latex');

%% 6c
Tmax = 3*2*pi;
figure
[xBDF2,tBDF2] = BDF2(fun, Jfun, x0, 5e-1, Tmax);
plot(tBDF2,xBDF2,'r','LineWidth',2);
hold on
plot(0:1e-2:Tmax,xorig(0:1e-2:Tmax),'b','LineWidth',2);
[xBDF2,tBDF2] = BDF2(fun, Jfun, x0, 1, Tmax);
hold on
plot(tBDF2,xBDF2,'m','LineWidth',2);
grid on
legend('h=3e-1','exact','h=1')
xlabel('t','Interpreter','latex');
ylabel('x(t)','Interpreter','latex');
title('BDF2 vs Exact Solution','Interpreter','latex');

%% 6d
H    = [0.2; 0.1; 0.05; 0.02; 0.01; 0.005; 0.001];
Tmax = 1;
for it = 1:length(H)
    h = H(it);
    [xAB2,tAB2] = AB2(fun, Jfun, x0, h, Tmax);
    errorAB2(it) = abs(xAB2(end)-xorig(1));
    [xBDF2,tBDF2] = BDF2(fun, Jfun, x0, h, Tmax);
    errorBDF2(it) = abs(xBDF2(end)-xorig(1));

end
figure
loglog(H,errorAB2,'r','LineWidth',2);
hold on
loglog(H,errorBDF2,'m','LineWidth',2);
grid on
legend('AB2','BDF2')
xlabel('Step Size(h)','Interpreter','latex');
ylabel('|error(1)|','Interpreter','latex');
title('Error(1) vs Step Size(h)','Interpreter','latex');