function [xout,tout] = AB2(fun, Jfun, x0, h, Tmax)

%% Settings
alpha = [1 -1 0];
beta  = [0 1.5 -0.5];

%% Solve ODE
[xout,tout] = solve_ode(alpha,beta,x0,Tmax,h,fun,Jfun);