function r=penalty_diffusive_coefficient_e(x,y,t)
%% 注意这里给的是不同于Beatrice书上的，在罚项也会带有扩散系数的情形，适用于Grote的波方程解法。
%r=sin(pi*x)*sin(pi*y)/10;
r=sin(x);
%r=sin(pi*x);
%r=exp(-x);
%r=exp(x*x);
%r=x;
%r=x^(1/2);
%r=x+1;
%r=x;
%r=1;