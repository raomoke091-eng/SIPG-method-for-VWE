function r=function_f(x,y,t)

%r=2*t*(1 - exp(x*(x - 1)))*(2*y - 1)^2*exp(t^2 + y^2 - y) + 4*t*(1 - exp(x*(x - 1)))*exp(t^2 + y^2 - y) + 2*t*(1 - exp(y*(y - 1)))*(2*x - 1)^2*exp(t^2 + x^2 - x) + 4*t*(1 - exp(y*(y - 1)))*exp(t^2 + x^2 - x) - (1 - exp(t^2))*(1 - exp(x*(x - 1)))*(2*y - 1)^2*exp(y*(y - 1)) - 2*(1 - exp(t^2))*(1 - exp(x*(x - 1)))*exp(y*(y - 1)) - (1 - exp(t^2))*(1 - exp(y*(y - 1)))*(2*x - 1)^2*exp(x*(x - 1)) - 2*(1 - exp(t^2))*(1 - exp(y*(y - 1)))*exp(x*(x - 1)) + 2*(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)))*(2*t^2 + 1)*exp(t^2);%a=1,b=1
%r=2*t*x*(1 - exp(x*(x - 1)))*(2*y - 1)^2*exp(t^2 + y^2 - y) + 4*t*x*(1 - exp(x*(x - 1)))*exp(t^2 + y^2 - y) + 2*t*x*(1 - exp(y*(y - 1)))*(2*x - 1)^2*exp(t^2 + x^2 - x) + 4*t*x*(1 - exp(y*(y - 1)))*exp(t^2 + x^2 - x) + 2*t*(1 - exp(y*(y - 1)))*(2*x - 1)*exp(t^2 + x^2 - x) - 2*x*(1 - exp(t^2))*(1 - exp(x*(x - 1)))*(2*y - 1)^2*exp(y*(y - 1)) - 4*x*(1 - exp(t^2))*(1 - exp(x*(x - 1)))*exp(y*(y - 1)) - 2*x*(1 - exp(t^2))*(1 - exp(y*(y - 1)))*(2*x - 1)^2*exp(x*(x - 1)) - 4*x*(1 - exp(t^2))*(1 - exp(y*(y - 1)))*exp(x*(x - 1)) - 2*(1 - exp(t^2))*(1 - exp(y*(y - 1)))*(2*x - 1)*exp(x*(x - 1)) + 2*(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)))*(2*t^2 + 1)*exp(t^2);%u同上一行 a=x,b=2x
%r=2*t*x*(1 - exp(x*(x - 1)))*(2*y - 1)^2*exp(t^2 + y^2 - y) + 4*t*x*(1 - exp(x*(x - 1)))*exp(t^2 + y^2 - y) + 2*t*x*(1 - exp(y*(y - 1)))*(2*x - 1)^2*exp(t^2 + x^2 - x) + 4*t*x*(1 - exp(y*(y - 1)))*exp(t^2 + x^2 - x) + 2*t*(1 - exp(y*(y - 1)))*(2*x - 1)*exp(t^2 + x^2 - x) - y*(1 - exp(t^2))*(1 - exp(x*(x - 1)))*(2*y - 1)^2*exp(y*(y - 1)) - 2*y*(1 - exp(t^2))*(1 - exp(x*(x - 1)))*exp(y*(y - 1)) - y*(1 - exp(t^2))*(1 - exp(y*(y - 1)))*(2*x - 1)^2*exp(x*(x - 1)) - 2*y*(1 - exp(t^2))*(1 - exp(y*(y - 1)))*exp(x*(x - 1)) - (1 - exp(t^2))*(1 - exp(x*(x - 1)))*(2*y - 1)*exp(y*(y - 1)) + 2*(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)))*(2*t^2 + 1)*exp(t^2); %u同上一行 a=x,b=y



%r=2*t*x*(1 - exp(x*(x - 1)))*(2*y - 1)^2*exp(t^2 + y^2 - y) + 4*t*x*(1 - exp(x*(x - 1)))*exp(t^2 + y^2 - y) + 2*t*x*(1 - exp(y*(y - 1)))*(2*x - 1)^2*exp(t^2 + x^2 - x) + 4*t*x*(1 - exp(y*(y - 1)))*exp(t^2 + x^2 - x) + 2*t*(1 - exp(y*(y - 1)))*(2*x - 1)*exp(t^2 + x^2 - x) - y*(1 - exp(t^2))*(1 - exp(x*(x - 1)))*(2*y - 1)^2*exp(y*(y - 1)) - 2*y*(1 - exp(t^2))*(1 - exp(x*(x - 1)))*exp(y*(y - 1)) - y*(1 - exp(t^2))*(1 - exp(y*(y - 1)))*(2*x - 1)^2*exp(x*(x - 1)) - 2*y*(1 - exp(t^2))*(1 - exp(y*(y - 1)))*exp(x*(x - 1)) - (1 - exp(t^2))*(1 - exp(x*(x - 1)))*(2*y - 1)*exp(y*(y - 1)) + 2*(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)))*(2*t^2 + 1)*exp(t^2);


%r=2*t*(pi^2*t^2 + 3*pi^2*t + 3)*sin(pi*x)*sin(pi*y);
%r=t*(2*pi^2*t^2*y*sin(pi*x)*sin(pi*y) - pi*t^2*sin(pi*x)*cos(pi*y) + 6*pi^2*t*x*sin(pi*x)*sin(pi*y) - 3*pi*t*sin(pi*y)*cos(pi*x) + 6*sin(pi*x)*sin(pi*y));%u同上一行 a=x,b=y


%r=(1.0*pi^2*sinh(0.5*t) + 0.25*cosh(0.5*t) + 2.0*pi^2*cosh(0.5*t))*sin(pi*x)*sin(pi*y);
%r=2*(2 - pi^2)*exp(-2*t)*sin(pi*x)*sin(pi*y);
%r=2*pi^2*(t + 1)*sin(pi*x)*sin(pi*y);
%r=(0.25 + 1.0*pi^2)*exp(-0.5*t)*sin(pi*x)*sin(pi*y);

%r=exp(-t)*sin(pi*x)*sin(pi*y);%a=1,b=1
%r=(4*sqrt(2)*pi^2*sin(pi*x)^2*sin(pi*y)*cos(pi*(y + 1/4)) + pi^2*sin(pi*x)^2 + 10*sin(pi*x)*sin(pi*y) - sqrt(2)*pi^2*sin(pi*y)*cos(pi*(y + 1/4)))*exp(-t)/10;%u同上一行 a=sin(pix)sin(piy)/10,b=sin(pix)cos(piy)/10
%r=(2*pi^2*sin(4*x)*sin(4*y)*sin(pi*x)*sin(pi*y) - 4*pi*sin(4*x)*sin(pi*x)*cos(4*y)*cos(pi*y) - 4*pi*sin(4*y)*sin(pi*y)*cos(4*x)*cos(pi*x) - 2*pi^2*sin(pi*x)^2*sin(pi*y)^2 + pi^2*sin(pi*x)^2*cos(pi*y)^2 + 10*sin(pi*x)*sin(pi*y) + pi^2*sin(pi*y)^2*cos(pi*x)^2)*exp(-t)/10;%u同上一行 a=sin(x),b=sin(y)
%r=(-2*pi^2*sin(x)*sin(pi*x)*sin(pi*y) + 2*pi^2*sin(y)*sin(pi*x)*sin(pi*y) + sin(pi*x)*sin(pi*y) - pi*sin(pi*x)*cos(y)*cos(pi*y) + pi*sin(pi*y)*cos(x)*cos(pi*x))*exp(-t);%u同上一行 a=sin(x),b=sin(y)   可以候选
%r=(-3*pi^2*sin(pi*x)^2*sin(pi*y) + 3*pi^2*sin(pi*x)*sin(pi*y)^2 + sin(pi*x)*sin(pi*y) - pi^2*sin(pi*x) + pi^2*sin(pi*y))*exp(-t);%u同上一行 a=sin(pix),b=sin(piy)   可以候选  和上一行效果差不多
%r=(-pi*(2*pi*sin(pi*x) + cos(pi*x))*exp(2*t + y)*sin(pi*y) + pi*(2*pi*sin(pi*y) + cos(pi*y))*exp(2*t + x)*sin(pi*x) + exp(2*t + x + y)*sin(pi*x)*sin(pi*y))*exp(-3*t - x - y);%u同上一行 a=exp(-x),b=exp(-y) 也可以候选 但效果相比上两行差一点
%r=(2*pi*x*exp(x^2)*sin(pi*y)*cos(pi*x) - 2*pi*y*exp(y^2)*sin(pi*x)*cos(pi*y) - 2*pi^2*exp(x^2)*sin(pi*x)*sin(pi*y) + 2*pi^2*exp(y^2)*sin(pi*x)*sin(pi*y) + sin(pi*x)*sin(pi*y))*exp(-t);%u同上一行 a=exp(x^2),b=exp(y^2)  目前最好 h=0.125的时候体现出sup-optimal
%r=(-2*pi^2*x*sin(pi*x)*sin(pi*y) + 2*pi^2*y*sin(pi*x)*sin(pi*y) + sin(pi*x)*sin(pi*y) - pi*sin(pi*(x - y)))*exp(-t);  %u同上一行 a=x b=y 目前最好
%r=(2*sqrt(x)*sqrt(y)*(-2*pi^2*sqrt(x) + 2*pi^2*sqrt(y) + 1)*sin(pi*x)*sin(pi*y) - pi*sqrt(x)*sin(pi*x)*cos(pi*y) + pi*sqrt(y)*sin(pi*y)*cos(pi*x))*exp(-t)/(2*sqrt(x)*sqrt(y));%u同上一行 a=sqrt(x),b=sqrt(y)
%r=(2*sqrt(x)*(-2*pi^2*sqrt(x)*sin(pi*x) + 2*pi^2*x*sin(pi*x) + sin(pi*x) -pi*cos(pi*x)) + pi*cos(pi*x))*exp(-t)*sin(pi*y)/(2*sqrt(x));%u同上一行 a=sqrt(x),b=x  效果不好
%r=(1 - 2*pi^2)*exp(-t)*sin(pi*x)*sin(pi*y);%u同上一行 a=x+1,b=x  效果差
%r=(2*pi^2*x*sin(pi*x) + sin(pi*x) - pi*cos(pi*x))*exp(-t)*sin(pi*y);%u同上一行  a=x,b=2x
%r=(-2*pi^2*x^2*sin(pi*x)*sin(pi*y) + 2*pi*x*sin(pi*y)*cos(pi*x) + 2*pi^2*y^2*sin(pi*x)*sin(pi*y) - 2*pi*y*sin(pi*x)*cos(pi*y) + sin(pi*x)*sin(pi*y))*exp(-t);%u同上一行 a=x,b=sin(y)
%r=(-2*pi*y*exp(y^2)*sin(pi*x)*cos(pi*y) + 2*pi^2*exp(y^2)*sin(pi*x)*sin(pi*y) - 2*pi^2*sin(x)*sin(pi*x)*sin(pi*y) + sin(pi*x)*sin(pi*y) + pi*sin(pi*y)*cos(x)*cos(pi*x))*exp(-t);%a=sinx,b=exp(y*y)



%r=(4*pi^2*t*exp(t^2) - 2*pi^2*(1 - exp(t^2)) + (4*t^2 + 2)*exp(t^2))*sin(pi*x)*sin(pi*y);
%r=(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1)))*exp(-t);
%r=-pi*(2*x - 1)^2*(exp(y*(y - 1)) - 1)*exp(x*(x - 1))*sinh(pi*t) - (2*x - 1)^2*(exp(y*(y - 1)) - 1)*exp(x*(x - 1))*cosh(pi*t) - pi*(2*y - 1)^2*(exp(x*(x - 1)) - 1)*exp(y*(y - 1))*sinh(pi*t) - (2*y - 1)^2*(exp(x*(x - 1)) - 1)*exp(y*(y - 1))*cosh(pi*t) + pi^2*(exp(x*(x - 1)) - 1)*(exp(y*(y - 1)) - 1)*cosh(pi*t) - 2*pi*(exp(x*(x - 1)) - 1)*exp(y*(y - 1))*sinh(pi*t) - 2*(exp(x*(x - 1)) - 1)*exp(y*(y - 1))*cosh(pi*t) - 2*pi*(exp(y*(y - 1)) - 1)*exp(x*(x - 1))*sinh(pi*t) - 2*(exp(y*(y - 1)) - 1)*exp(x*(x - 1))*cosh(pi*t);
%r=-0.1*sin(t)*sin(pi*(10*x + 10*y)) + 20.0*pi^2*sin(t)*sin(pi*(10*x + 10*y)) - 0.125*pi^2*sin(t)*cos(pi*(0.25*x + 0.25*y)) + 20.0*pi^2*sin(pi*(10*x + 10*y))*cos(t) - cos(t)*cos(pi*(0.25*x + 0.25*y)) + 0.125*pi^2*cos(t)*cos(pi*(0.25*x + 0.25*y));
%r=(-sin(t) + 200*sqrt(2)*pi^2*sin(t + pi/4))*sin(pi*(10*x + 10*y));
%r=2*(2 - pi^2)*exp(-2*t)*sin(pi*(x + y));
%r=(pi^2*(1 - exp(x*(x - 1))) + (2*x - 1)^2*exp(x*(x - 1)) + 6*exp(x*(x - 1)) - 4)*exp(-2*t)*sin(pi*y);%a=1,b=1
%r=4*(exp(x*(x - 1)) - 1)*exp(-2*t)*sin(pi*y);% u倒数第二行 a=x,b=2x



%r=(4*(1 - exp(x*(x - 1)))*(1 - exp(y*(y - 1))) - (1 - exp(x*(x - 1)))*(2*y - 1)^2*exp(y*(y - 1)) - 2*(1 - exp(x*(x - 1)))*exp(y*(y - 1)) - (1 - exp(y*(y - 1)))*(2*x - 1)^2*exp(x*(x - 1)) - 2*(1 - exp(y*(y - 1)))*exp(x*(x - 1)))*exp(-2*t);



%r=exp(-t)*sin(2*pi*x)*sin(2*pi*y);


%r=(2*pi^2*x*sin(pi*x) + sin(pi*x) - pi*cos(pi*x))*exp(-t)*sin(pi*y);
r=(-2*pi*y*exp(y^2)*sin(pi*x)*cos(pi*y) + 2*pi^2*exp(y^2)*sin(pi*x)*sin(pi*y) - 2*pi^2*sin(x)*sin(pi*x)*sin(pi*y) + sin(pi*x)*sin(pi*y) + pi*sin(pi*y)*cos(x)*cos(pi*x))*exp(-t);
