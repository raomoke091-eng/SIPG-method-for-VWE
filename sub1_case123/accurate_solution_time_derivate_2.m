function r=accurate_solution_time_derivate_2(x,y,t)

%r=sin(10*x)*sin(10*y)*exp(-t);
r=2*(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)));
%r=sin(5*pi*x)*sin(5*pi*y)*cosh(t);
%r=exp(-t)*sin(2*pi*x)*sin(2*pi*y);