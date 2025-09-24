function r=accurate_solution(x,y,t)

%r=sin(pi*x)*sin(pi*y)*exp(-t);
%r=sin(pi*x)*sin(pi*y)*t^2;
r=sin(pi*x)*sin(pi*y)*t^3;
%r=sin(pi*x)*sin(pi*y)*sin(t);
%r=sin(pi*x)*sin(pi*y)*tanh(t);
%r=sin(pi*x)*sin(pi*y)*t^4;
%r=sin(pi*x)*sin(pi*y)*t;
%r=sin(pi*x)*sin(pi*y)*exp(-3*t);
%r=sin(pi*x)*sin(pi*y)*log(t);
%r=cos(t)*cos(0.25*pi*(x+y))+0.1*sin(t)*sin(10*pi*(x+y));