function [xout,tout] = BDF2(fun, Jfun, x0, h, Tmax)

%% Settings
alpha = [1 -1.3333 0.3333];
beta  = [0.6667 0 0];

%% Solve ODE
[xout,tout] = solve_ode(alpha,beta,x0,Tmax,h,fun,Jfun);