function [x,t,counter,dx] = solve_ode(alpha, beta, x0, Tmax, h, fun, varargin)
%
if nargin > 6
    Jf = varargin{1};
end

%% Settings
eps = 1e-8;
max_iter = 1e1;

% Initializations
alpha = alpha./alpha(1); % normalize alpha
alpha = reshape(alpha,1,length(alpha));
beta  = reshape(beta,1,length(beta));
k    = length(alpha)-1;
t = 0:h:Tmax;

num_steps = length(t);
n = length(x0);
x = zeros(n,num_steps);
x(:,1)= x0;
% x(:,2)= cos(t(2));

%% Initialize with Backward Euler
if k > 1
    x(:,1:k) = backward_euler(fun,Jf,x0,h, h*(k-1));
end

if beta(1)==0 % If Explicit Method
    for i = k+1:num_steps
        u      = t(i);
        x(:,i) = -x(:,i-1:-1:i-k)*alpha(2:end)' + h * fun(x(:,i-1:-1:i-length(beta)+1),u)*beta(2:end)';
    end
else
    % Implicit Method
    % Construct Jacobian and RHS
    for i = k+1:num_steps
        u = t(i);
        b = x(:,i-1:-1:i-k)*alpha(2:end)' - h * fun(x(:,i-1:-1:i-length(beta)+1),u)*beta(2:end)';
        Jacobianfun = @(x) eye(n) - h * beta(1) * Jf(x);
        RHSfun      = @(x) (x - h * beta(1) * fun(x,u) + b);
        [x(:,i),counter,dx] = solvenewton(RHSfun, Jacobianfun, x(:,i-1), max_iter, eps);
    end
end