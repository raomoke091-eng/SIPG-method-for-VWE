function r=accurate_solution_x_derivative_2(x,y,t)

%r=-100*sin(10*x)*sin(10*y)*exp(-t);
r=t^2*((2*x - 1)^2 + 2)*(exp(y*(y - 1)) - 1)*exp(x*(x - 1));
%r=-25*pi^2*sin(5*pi*x)*sin(5*pi*y)*cosh(t);
%r=-4*pi^2*exp(-t)*sin(2*pi*x)*sin(2*pi*y);