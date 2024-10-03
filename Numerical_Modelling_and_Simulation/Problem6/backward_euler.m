function [xout,tout] = backward_euler(fun, Jfun, x0, h, Tmax)

%% Settings
alpha = [1 -1];
beta  = [1  0];

%% Solve ODE
[xout,tout] = solve_ode(alpha, beta, x0, Tmax, h, fun, Jfun);