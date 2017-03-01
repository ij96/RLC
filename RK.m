function [t,x,y] = RK(method,dx,dy,xi,yi,ti,tf,h)
    % Runge-Kutta method for solving the 2nd order ODE y'' = f(t,x,y)
    % This is achieved by writing the 2nd order ODE as two coupled 1st order ODEs:
    %   dx = y = x'
    %   dy = y' = x''
    % then solving both together using a specified RK method.
    %
    % Parameters:
    %   method: name of the RK method to be used (case-insensitive)
    %   dx, dy: coupled 1st order ODEs equivalent of the 2nd order ODE y'' = f(t,x,y)
    %       dx = y = x'
    %       dy = y' = x''
    %   xi: initial value of x
    %   yi: initial value of y
    %   ti: initial value of t
    %   tf: final value of t
    %   h: time interval size
    
    % Generate the necessary matrices:
    %   s: number of stages = order of the method
    %   a: Runge-Kutta matrix
    %   b: weights
    %   c: nodes
    if(strcmpi(method,'euler'))          % Euler's method, RK1
        s = 1;
        a = 0;
        b = 1;
        c = 0;
    elseif(strcmpi(method,'midpoint'))   % Explicit midpoint method, RK2
        s = 2;
        a = [0    0;
             1/2  0];
        b = [0    1];
        c = [0    1/2];
    elseif(strcmpi(method,'heun'))       % Heun's method, RK2
        s = 2;
        a = [0    0;
             1    0];
        b = [1/2  1/2];
        c = [0    1];
    elseif(strcmpi(method,'ralston'))    % Ralston's method, RK2
        s = 2;
        a = [0    0;
             2/3    0];
        b = [1/4  3/4];
        c = [0    2/3];
    elseif(strcmpi(method,'rk3'))        % Kutta's 3rd order method, RK3
        s = 3;
        a = [0    0    0;
             1/2  0    0;
             -1   2    0];
        b = [1/6  2/3  1/6];
        c = [0    1/2  1];
    elseif(strcmpi(method,'rk4'))        % Classic 4th order method, RK4
        s = 4;
        a = [0    0    0    0;
             1/2  0    0    0;
             0    1/2  0    0;
             0    0    1    0];
        b = [1/6  1/3  1/3  1/6];
        c = [0    1/2  1/2  1];
    elseif(strcmpi(method,'3/8'))        % 3/8 4th order method, RK4
        s = 4;
        a = [0    0    0    0;
             1/3  0    0    0;
             -1/3 1    0    0;
             1    -1   1    0];
        b = [1/8  3/8  3/8  3/8];
        c = [0    1/3  2/3  1];
    else
        error('Method cannot be found.');
    end
    
    % Applying numerical method
    N = round((tf-ti)/h);       % total number of steps
    t = zeros(1,N);             % [1 x N] matrix of 0's
    x = zeros(1,N);
    y = zeros(1,N);

    % Set initial conditions
    t(1) = ti;
    x(1) = xi;
    y(1) = yi;
    
    % Iteration
    for n = 1:N-1
        k = zeros(1,s);
        l = zeros(1,s);
        for j = 1:s
            tt = t(n) + h * c(j);
            xx = x(n) + dot(a(j,:),k);
            yy = y(n) + dot(a(j,:),l);
            k(j) = h * dx(tt, xx, yy);
            l(j) = h * dy(tt, xx, yy);
        end
        t(n+1) = t(n) + h;
        x(n+1) = x(n) + dot(b,k);
        y(n+1) = y(n) + dot(b,l);
    end