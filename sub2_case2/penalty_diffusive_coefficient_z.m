function r=penalty_diffusive_coefficient_z(x,y,t)
%% 注意这里给的是不同于Beatrice书上的，在罚项也会带有扩散系数的情形，适用于Grote的波方程解法。
%r=exp(x);
%r=x*x;
%r=1;
%r=sin(pi*x)*cos(pi*y)/10;
%r=sin(4*x)*sin(4*y)/10;
%r=sin(y);
%r=sin(pi*y);
%r=exp(-y);
r=exp(y*y);
%r=y;
%r=y^(1/2);
%r=2*x;
%r=y;
%r=1;