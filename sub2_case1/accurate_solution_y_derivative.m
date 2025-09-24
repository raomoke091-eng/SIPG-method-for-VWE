function r=accurate_solution_y_derivative(x,y,t)

%r=(1 - exp(t^2))*(1 - exp(x*(x - 1)))*(2*y - 1)*exp(y*(y - 1));
%r=pi*t^3*sin(pi*x)*cos(pi*y);
%r=pi*sin(pi*x)*cos(pi*y)*cosh(0.5*t);
%r=pi*exp(-2*t)*sin(pi*x)*cos(pi*y);
%r=pi*t*sin(pi*x)*cos(pi*y);
%r=pi*exp(-0.5*t)*sin(pi*x)*cos(pi*y);
%r=pi*exp(-t)*sin(pi*x)*cos(pi*y);
%r=pi*(exp(t^2) - 1)*sin(pi*x)*cos(pi*y);
%r=(2*y - 1)*(exp(x*(x - 1)) - 1)*exp(-t + y^2 - y);
%r=(2*y - 1)*(exp(x*(x - 1)) - 1)*exp(y*(y - 1))*cosh(pi*t);
%r=pi*(1.0*sin(t)*cos(pi*(10*x + 10*y)) - 0.25*sin(pi*(0.25*x + 0.25*y))*cos(t));
%r=10*pi*sin(t)*cos(pi*(10*x + 10*y));
%r=pi*exp(-2*t)*cos(pi*(x + y));
%r=pi*(exp(x*(x - 1)) - 1)*exp(-2*t)*cos(pi*y);
r=(2*y - 1)*(exp(x*(x - 1)) - 1)*exp(-2*t + y^2 - y);