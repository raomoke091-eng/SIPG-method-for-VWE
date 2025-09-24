function r=accurate_solution_x_derivative_y_derivative(x,y,t)

%r=(2*x - 1)*(2*y - 1)*(exp(t^2) - 1)*exp(x^2 - x + y^2 - y);
%r=pi^2*t^3*cos(pi*x)*cos(pi*y);
%r=pi^2*cos(pi*x)*cos(pi*y)*cosh(0.5*t);
%r=pi^2*exp(-2*t)*cos(pi*x)*cos(pi*y);
%r=pi^2*t*cos(pi*x)*cos(pi*y);
%r=pi^2*exp(-0.5*t)*cos(pi*x)*cos(pi*y);
r=pi^2*exp(-t)*cos(pi*x)*cos(pi*y);
%r=pi^2*(exp(t^2) - 1)*cos(pi*x)*cos(pi*y);
%r=(2*x - 1)*(2*y - 1)*exp(-t + x^2 - x + y^2 - y);
%r=(2*x - 1)*(2*y - 1)*exp(x^2 - x + y^2 - y)*cosh(pi*t);
%r=pi^2*(-10.0*sin(t)*sin(pi*(10*x + 10*y)) - 0.0625*cos(t)*cos(pi*(0.25*x + 0.25*y)));
%r=-100*pi^2*sin(t)*sin(pi*(10*x + 10*y));
%r=-pi^2*exp(-2*t)*sin(pi*(x + y));
%r=pi*(2*x - 1)*exp(-2*t + x^2 - x)*cos(pi*y);
%r=(2*x - 1)*(2*y - 1)*exp(-2*t + x^2 - x + y^2 - y);
%r=4*pi^2*exp(-t)*cos(2*pi*x)*cos(2*pi*y);