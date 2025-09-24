function r=accurate_solution(x,y,t)

%r=(exp(x*x-x)-1)*(exp(y*y-y)-1)*(exp(t*t)-1);
%r=sin(pi*x)*sin(pi*y)*t^3;
%r=sin(pi*x)*sin(pi*y)*cosh(0.5*t);%可以候选
%r=sin(pi*x)*sin(pi*y)*exp(-2*t);%可以候选
%r=sin(pi*x)*sin(pi*y)*t;
%r=sin(pi*x)*sin(pi*y)*exp(-0.5*t);%可以候选
r=sin(pi*x)*sin(pi*y)*exp(-t);
%r=sin(pi*x)*sin(pi*y)*(exp(t^2)-1);
%r=(exp(x*x-x)-1)*(exp(y*y-y)-1)*exp(-t);
%r=(exp(x*x-x)-1)*(exp(y*y-y)-1)*cosh(0.5*t);
%r=cos(t)*cos(0.25*pi*(x+y))+0.1*sin(t)*sin(10*pi*(x+y));
%r=sin(t)*sin(10*pi*(x+y));
%r=sin(pi*(x+y))*exp(-2*t);
%r=(exp(x*x-x)-1)*sin(pi*y)*exp(-2*t);%可以候选 目前最好
%r=(exp(x*x-x)-1)*(exp(y*y-y)-1)*exp(-2*t);%可以候选 次好
%r=sin(2*pi*x)*sin(2*pi*y)*exp(-t);