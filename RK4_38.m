function [x,y] = RK4_38(dx,dy,ti,xi,yi,h)
    % Rnge-Kutta 3/8 method for solving the 2st order ODE y'' = f(t,x,y)
    % This function gives the next values (x,y) from the previous values (ti,xi,yi)
    %
    % dx = x' = y, dy = y' = x''
    k1 = h * dx(ti,       xi,     	   yi);
    l1 = h * dy(ti,       xi,          yi);
    k2 = h * dx(ti+h/3,   xi+k1/3,     yi+l1/3);
    l2 = h * dy(ti+h/3,   xi+k1/3,     yi+l1/3);
    k3 = h * dx(ti+2*h/3, xi-k1/3+k2,  yi-l1/3+l2);
    l3 = h * dy(ti+2*h/3, xi-k1/3+k2,  yi-l1/3+l2);
    k4 = h * dx(ti+h,     xi+k1-k2+k3, yi+l1-l2+l3);
    l4 = h * dy(ti+h,     xi+k1-k2+k3, yi+l1-l2+l3);
    x  = xi + 1 / 8 * (k1 + 3*k2 + 3*k3 + k4);
    y  = yi + 1 / 8 * (l1 + 3*l2 + 3*l3 + l4);