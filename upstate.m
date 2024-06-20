function xnext = upstate(xdot,xi,step)

% xnext = euler2(xdot,x,h)
%
% xdot  - dx/dt(k) = f(x(k)
% x     - x(k)
% xnext - x(k+1)
% h     - step size
xnext = xi + step*xdot;

