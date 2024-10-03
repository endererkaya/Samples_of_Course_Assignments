function [xval,counter,dx] = solvenewton(fun, Jfun, xinit, max_iter, eps)

% Settings
% eps      = 1e-8;

% Initializations
counter  = 0;
xval     = xinit;
dx = [];
% Newton Algorithm
while (counter < max_iter)
    counter = counter + 1;

    % Newton Update
    funval  = fun(xval); % evaluate F(x)
    JFt     = Jfun(xval); % evaluate JacobianF(x)
    deltax  = -JFt \ funval;
    xval    = xval + deltax;
    
%     xxplot  = [xxplot xval];
    dx      = [dx deltax];
    % Convergence Check
    if (norm(deltax) < (eps * norm(xval))) && (norm(funval) < eps)
        break;
    end
end