function r=accurate_solution_x_derivative(x,y,t)

%r=(1 - exp(t^2))*(1 - exp(y*(y - 1)))*(2*x - 1)*exp(x*(x - 1));
%r=pi*t^3*sin(pi*y)*cos(pi*x);
%r=pi*sin(pi*y)*cos(pi*x)*cosh(0.5*t);
%r=pi*exp(-2*t)*sin(pi*y)*cos(pi*x);
%r=pi*t*sin(pi*y)*cos(pi*x);
%r=pi*exp(-0.5*t)*sin(pi*y)*cos(pi*x);
r=pi*exp(-t)*sin(pi*y)*cos(pi*x);
%r=pi*(exp(t^2) - 1)*sin(pi*y)*cos(pi*x);
%r=(2*x - 1)*(exp(y*(y - 1)) - 1)*exp(-t + x^2 - x);
%r=(2*x - 1)*(exp(y*(y - 1)) - 1)*exp(x*(x - 1))*cosh(pi*t);
%r=pi*(1.0*sin(t)*cos(pi*(10*x + 10*y)) - 0.25*sin(pi*(0.25*x + 0.25*y))*cos(t));
%r=10*pi*sin(t)*cos(pi*(10*x + 10*y));
%r=pi*exp(-2*t)*cos(pi*(x + y));
%r=(2*x - 1)*exp(-2*t + x^2 - x)*sin(pi*y);
%r=(2*x - 1)*(exp(y*(y - 1)) - 1)*exp(-2*t + x^2 - x);
%r=2*pi*exp(-t)*sin(2*pi*y)*cos(2*pi*x);