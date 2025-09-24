function r=accurate_solution_y_derivative(x,y,t)

%r=10*sin(10*x)*cos(10*y)*exp(-t);
r=t^2*(2*y - 1)*(exp(x*(x - 1)) - 1)*exp(y*(y - 1));
%r=5*pi*sin(5*pi*x)*cos(5*pi*y)*cosh(t);
%r=2*pi*exp(-t)*sin(2*pi*x)*cos(2*pi*y);