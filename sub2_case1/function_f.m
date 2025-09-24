function r=function_f(x,y,t)

%r=2*t*(1 - exp(x*(x - 1)))*(2*y - 1)^2*exp(t^2 + y^2 - y) + 4*t*(1 - exp(x*(x - 1)))*exp(t^2 + y^2 - y) + 2*t*(1 - exp(y*(y - 1)))*(2*x - 1)^2*exp(t^2 + x^2 - x) + 4*t*(1 - exp(y*(y - 1)))*exp(t^2 + x^2 - x) - (1 - exp(t^2))*(1 - exp(x*(x - 1)))*(2*y - 1)^2*exp(y*(y - 1)) - 2*(1 - exp(t^2))*(1 - exp(x*(x - 1)))*exp(y*(y - 1)) - (1 - exp(t^2))*(1 - exp(y*(y - 1)))*(2*x - 1)^2*exp(x*(x - 1)) - 2*(1 - exp(t^2))*(1 - exp(y*(y - 1)))*exp(x*(x - 1)) + 2*(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)))*(2*t^2 + 1)*exp(t^2);
%r=2*t*(pi^2*t^2 + 3*pi^2*t + 3)*sin(pi*x)*sin(pi*y);
%r=(1.0*pi^2*sinh(0.5*t) + 0.25*cosh(0.5*t) + 2.0*pi^2*cosh(0.5*t))*sin(pi*x)*sin(pi*y);
%r=2*(2 - pi^2)*exp(-2*t)*sin(pi*x)*sin(pi*y);
%r=2*pi^2*(t + 1)*sin(pi*x)*sin(pi*y);
%r=(0.25 + 1.0*pi^2)*exp(-0.5*t)*sin(pi*x)*sin(pi*y);
%r=exp(-t)*sin(pi*x)*sin(pi*y);
%r=(4*pi^2*t*exp(t^2) - 2*pi^2*(1 - exp(t^2)) + (4*t^2 + 2)*exp(t^2))*sin(pi*x)*sin(pi*y);
%r=(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)))*exp(-t);
%r=-pi*(2*x - 1)^2*(exp(y*(y - 1)) - 1)*exp(x*(x - 1))*sinh(pi*t) - (2*x - 1)^2*(exp(y*(y - 1)) - 1)*exp(x*(x - 1))*cosh(pi*t) - pi*(2*y - 1)^2*(exp(x*(x - 1)) - 1)*exp(y*(y - 1))*sinh(pi*t) - (2*y - 1)^2*(exp(x*(x - 1)) - 1)*exp(y*(y - 1))*cosh(pi*t) + pi^2*(exp(x*(x - 1)) - 1)*(exp(y*(y - 1)) - 1)*cosh(pi*t) - 2*pi*(exp(x*(x - 1)) - 1)*exp(y*(y - 1))*sinh(pi*t) - 2*(exp(x*(x - 1)) - 1)*exp(y*(y - 1))*cosh(pi*t) - 2*pi*(exp(y*(y - 1)) - 1)*exp(x*(x - 1))*sinh(pi*t) - 2*(exp(y*(y - 1)) - 1)*exp(x*(x - 1))*cosh(pi*t);
%r=-0.1*sin(t)*sin(pi*(10*x + 10*y)) + 20.0*pi^2*sin(t)*sin(pi*(10*x + 10*y)) - 0.125*pi^2*sin(t)*cos(pi*(0.25*x + 0.25*y)) + 20.0*pi^2*sin(pi*(10*x + 10*y))*cos(t) - cos(t)*cos(pi*(0.25*x + 0.25*y)) + 0.125*pi^2*cos(t)*cos(pi*(0.25*x + 0.25*y));
%r=(-sin(t) + 200*sqrt(2)*pi^2*sin(t + pi/4))*sin(pi*(10*x + 10*y));
%r=2*(2 - pi^2)*exp(-2*t)*sin(pi*(x + y));
%r=(pi^2*(1 - exp(x*(x - 1))) + (2*x - 1)^2*exp(x*(x - 1)) + 6*exp(x*(x - 1)) - 4)*exp(-2*t)*sin(pi*y);
r=(4*(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1))) - (1 - exp(x*(x - 1)))*(2*y - 1)^2*exp(y*(y - 1)) - 2*(1 - exp(x*(x - 1)))*exp(y*(y - 1)) - (1 - exp(y*(y - 1)))*(2*x - 1)^2*exp(x*(x - 1)) - 2*(1 - exp(y*(y - 1)))*exp(x*(x - 1)))*exp(-2*t);