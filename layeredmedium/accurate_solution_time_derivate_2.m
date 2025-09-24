function r=accurate_solution_time_derivate_2(x,y,t)

%r=sin(pi*x)*sin(pi*y)*exp(-t);
%r=sin(pi*x)*sin(pi*y)*2;
r=sin(pi*x)*sin(pi*y)*t*6;
%r=-sin(pi*x)*sin(pi*y)*sin(t);
%r=-2*sin(pi*x)*sin(pi*y)*(sech(t))^2*tanh(t);
%r=sin(pi*x)*sin(pi*y)*t^2*12;
%r=0;
%r=sin(pi*x)*sin(pi*y)*exp(-3*t)*9;
%r=sin(pi*x)*sin(pi*y)*(-1/(t^2));
%r=-0.1 * sin(t) * sin(10 * pi * (x + y)) - cos(t) * cos(0.25 * pi * (x + y));