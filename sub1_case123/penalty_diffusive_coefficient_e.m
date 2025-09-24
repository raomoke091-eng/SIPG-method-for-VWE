function r=penalty_diffusive_coefficient_e(x,y,t)
%% 注意这里给的是不同于Beatrice书上的，在罚项也会带有扩散系数的情形，适用于Grote的波方程解法。
%r=1/20;
%r=0.00102;
%r=sin(pi*x)*sin(pi*y)/8;
%r=0.0001;
%r=0.00093;
r=1;

%r=0.056;  %干砂岩阻尼系数
%r=0.2;  %饱和水砂岩阻尼系数
%r=0.0112;  %气体饱和层

%suboptimal test
%r=1;
%r=sin(x+y);