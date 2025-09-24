function r=accurate_solution_time_derivate(x,y,t)

%r=-sin(pi*x)*sin(pi*y)*exp(-t);
%r=sin(pi*x)*sin(pi*y)*t*2;
r=sin(pi*x)*sin(pi*y)*t^2*3;
%r=sin(pi*x)*sin(pi*y)*cos(t);
%r=sin(pi*x)*sin(pi*y)*(sech(t))^2;
%r=sin(pi*x)*sin(pi*y)*t^3*4;
%r=sin(pi*x)*sin(pi*y);
%r=sin(pi*x)*sin(pi*y)*exp(-3*t)*(-3);
%r=sin(pi*x)*sin(pi*y)*(1/t);
%r=-sin(t) * cos(0.25 * pi * (x + y)) + 0.1 * sin(10 * pi * (x + y)) * cos(t);
