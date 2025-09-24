function r=accurate_solution_time_derivate_2(x,y,t)

%r=2*(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)))*(2*t^2 + 1)*exp(t^2);
%r=6*t*sin(pi*x)*sin(pi*y);
%r=0.25*sin(pi*x)*sin(pi*y)*cosh(0.5*t);
%r=4*exp(-2*t)*sin(pi*x)*sin(pi*y);
%r=0;
%r=0.25*exp(-0.5*t)*sin(pi*x)*sin(pi*y);
r=exp(-t)*sin(pi*x)*sin(pi*y);
%r=2*(2*t^2 + 1)*exp(t^2)*sin(pi*x)*sin(pi*y);
%r=(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)))*exp(-t);
%r=pi^2*(exp(x*(x - 1)) - 1)*(exp(y*(y - 1)) - 1)*cosh(pi*t);
%r=0.1*sin(t)*sin(pi*(10*x + 10*y)) - cos(t)*cos(pi*(0.25*x + 0.25*y));
%r=-sin(t)*sin(pi*(10*x + 10*y));
%r=4*exp(-2*t)*sin(pi*(x + y));
%r=4*(exp(x*(x - 1)) - 1)*exp(-2*t)*sin(pi*y);
%r=4*(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)))*exp(-2*t);
%r=exp(-t)*sin(2*pi*x)*sin(2*pi*y);