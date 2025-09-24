function r=accurate_solution_x_derivative(x,y,t)

%r=10*cos(10*x)*sin(10*y)*exp(-t);
r=t^2*(2*x - 1)*(exp(y*(y - 1)) - 1)*exp(x*(x - 1));
%r=5*pi*sin(5*pi*y)*cos(5*pi*x)*cosh(t);
%r=2*pi*exp(-t)*sin(2*pi*y)*cos(2*pi*x);