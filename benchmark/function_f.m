function r=function_f(x,y,t)
f_0 = 100; % 可以根据实际情况设置主频赫兹数
x_s = 0; %135 % x 中心点
y_s = 880; %800 % y 中心点
R = sqrt((x-x_s)^2+(y-y_s)^2);
t_0 = 1/f_0;
Iter_refinement=5;
h = (960-800)*2^(-Iter_refinement);

%时间传播函数
if t <= 2 * t_0
    F_t = -2 * pi^2 * f_0^2 * (t - t_0) * exp(-pi^2 * f_0^2 * (t - t_0)^2);
else
    F_t = 0;
end

%空间传播函数
a = 5*h;
if R <= a
    indicator_value = 1;
else
    indicator_value = 0;
end

g_r = (1 - (R/a)^2)* indicator_value* (y-y_s)/R;

%合成
r=F_t*g_r;

