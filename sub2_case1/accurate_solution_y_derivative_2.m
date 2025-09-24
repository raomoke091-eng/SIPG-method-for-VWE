function r=accurate_solution_y_derivative_2(x,y,t)

%r=(1 - exp(t^2))*(1 - exp(x*(x - 1)))*((2*y - 1)^2 + 2)*exp(y*(y - 1));
%r=-pi^2*t^3*sin(pi*x)*sin(pi*y);
%r=-pi^2*sin(pi*x)*sin(pi*y)*cosh(0.5*t);
%r=-pi^2*exp(-2*t)*sin(pi*x)*sin(pi*y);
%r=-pi^2*t*sin(pi*x)*sin(pi*y);
%r=-pi^2*exp(-0.5*t)*sin(pi*x)*sin(pi*y);
%r=-pi^2*exp(-t)*sin(pi*x)*sin(pi*y);
%r=pi^2*(1 - exp(t^2))*sin(pi*x)*sin(pi*y);
%r=((2*y - 1)^2 + 2)*(exp(x*(x - 1)) - 1)*exp(-t + y*(y - 1));
%r=((2*y - 1)^2 + 2)*(exp(x*(x - 1)) - 1)*exp(y*(y - 1))*cosh(pi*t);
%r=pi^2*(-10.0*sin(t)*sin(pi*(10*x + 10*y)) - 0.0625*cos(t)*cos(pi*(0.25*x + 0.25*y)));
%r=-100*pi^2*sin(t)*sin(pi*(10*x + 10*y));
%r=-pi^2*exp(-2*t)*sin(pi*(x + y));
%r=pi^2*(1 - exp(x*(x - 1)))*exp(-2*t)*sin(pi*y);
r=((2*y - 1)^2 + 2)*(exp(x*(x - 1)) - 1)*exp(-2*t + y*(y - 1));