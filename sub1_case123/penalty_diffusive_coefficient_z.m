function r=penalty_diffusive_coefficient_z(x,y,t)
%% 注意这里给的是不同于Beatrice书上的，在罚项也会带有扩散系数的情形，适用于Grote的波方程解法。

%r=1/20;
%r=3.06*10^(-9);
%r=sin(pi*x)*sin(pi*y)/8;
%r=4.2*10^(-10);
%r=5.12*10^(-8);
r=1;
%r=x;

%r=1190^2;  %干砂岩扩散系数
%r=1470^2;  %饱和水砂岩扩散系数
%r=630^2;  %气体饱和层扩散系数

%suboptimal test
%r=x;
%r=exp(y*y);
%r=cos(x+y);