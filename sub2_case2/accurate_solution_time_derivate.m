function r=accurate_solution_time_derivate(x,y,t)

%r=2*t*(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)))*exp(t^2);
%r=3*t^2*sin(pi*x)*sin(pi*y);
%r=0.5*sin(pi*x)*sin(pi*y)*sinh(0.5*t);
%r=-2*exp(-2*t)*sin(pi*x)*sin(pi*y);
%r=sin(pi*x)*sin(pi*y);
%r=-0.5*exp(-0.5*t)*sin(pi*x)*sin(pi*y);
r=-exp(-t)*sin(pi*x)*sin(pi*y);
%r=2*t*exp(t^2)*sin(pi*x)*sin(pi*y);
%r=-(exp(x*(x - 1)) - 1)*(exp(y*(y - 1)) - 1)*exp(-t);
%r=pi*(exp(x*(x - 1)) - 1)*(exp(y*(y - 1)) - 1)*sinh(pi*t);
%r=-sin(t)*cos(pi*(0.25*x + 0.25*y)) + 0.1*sin(pi*(10*x + 10*y))*cos(t);
%r=sin(pi*(10*x + 10*y))*cos(t);
%r=-2*exp(-2*t)*sin(pi*(x + y));
%r=2*(1 - exp(x*(x - 1)))*exp(-2*t)*sin(pi*y);
%r=-2*(exp(x*(x - 1)) - 1)*(exp(y*(y - 1)) - 1)*exp(-2*t);
%r=-exp(-t)*sin(2*pi*x)*sin(2*pi*y);