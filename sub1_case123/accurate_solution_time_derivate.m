function r=accurate_solution_time_derivate(x,y,t)

%r=-exp(-t)*sin(10*x)*sin(10*y);
r=2*t*(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)));
%r=sin(5*pi*x)*sin(5*pi*y)*sinh(t);
%r=-exp(-t)*sin(2*pi*x)*sin(2*pi*y);