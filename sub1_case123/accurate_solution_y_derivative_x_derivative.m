function r=accurate_solution_y_derivative_x_derivative(x,y,t)

%r=100*cos(10*x)*cos(10*y)*exp(-t);
r=t^2*(2*x - 1)*(2*y - 1)*exp(x^2 - x + y^2 - y);
%r=25*pi^2*cos(5*pi*x)*cos(5*pi*y)*cosh(t);
%r=4*pi^2*exp(-t)*cos(2*pi*x)*cos(2*pi*y);