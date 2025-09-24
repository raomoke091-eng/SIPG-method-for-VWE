function r=accurate_solution(x,y,t)

%r=sin(10*x)*sin(10*y)*exp(-t); %case1
r=(exp(x*x-x)-1)*(exp(y*y-y)-1)*t^2; %case2
%r=sin(5*pi*x)*sin(5*pi*y)*cosh(t); %case3
%r=sin(2*pi*x)*sin(2*pi*y)*exp(-t);
