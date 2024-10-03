%% ELEC 518 HW4 Problem 5
%
% Ender Erkaya
% 
% January 2022
%
%% Take Input
% k = input('Please enter k values 1,2,3,..  \n');
x = zeros(2*k+2,1);

% a = input('Enter array elements [ ] around them \n');
x(1:length(a)) = a; %array coefficients
alpha = x(1:k+1);
beta  = x(k+2:2*k+2);

% h = input('Step size?\n');

N = 1e3;
theta = linspace(0,2*pi,N);

zvals = exp(1i.*theta);

% Create Vandermonde Matrix
Zvan  = zeros(N,k+1);
for c = 1:k+1
    Zvan(:,c)=zvals.^(k-c+1);
end

% Form Equation
rhs = Zvan*alpha;
lhs = h*(Zvan*beta); % lambda coefficient
lambdas = rhs./lhs;

%% Create Visual
len = linspace(-2*max(abs(lambdas)),2*max(abs(lambdas)),500);
% r = 0:1e-2:2*max(abs(lambdas));
% theta = linspace(0,2*pi,1e3);
% [rr, tt] = meshgrid(r,theta);
% x = rr.*cos(tt);
% y = rr.*sin(tt);
[x, y] = meshgrid(len, len);

z = zeros(size(x));
for r = 1:size(x,1)
    r
    for s = 1:size(y,2)
        temp = x(r,s) + y(r,s)*1i;
        
        % Determine if stable
        if max(abs(roots(alpha-temp*h*beta)))<1
            z(r,s)= -0.2;
        else
            z(r,s)= 1;
        end
    end
end

contourf(x,y,z,'k','LineWidth',2);
axis equal
hold on
line([min(len) max(len)],[0 0],'Color','black','LineStyle','--');
hold on
line([0 0],[min(len) max(len)],'Color','black','LineStyle','--')
grid on
caxis([-1, 1]);
colormap('hot');
% title(strcat('Stability Region, h =',string(h)),'Interpreter','latex')
xlabel('Re($\lambda$)','Interpreter','latex');
ylabel('Im($\lambda$)','Interpreter','latex');