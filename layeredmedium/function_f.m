function r=function_f(x,y,t)
f_0 = 10.2; % 可以根据实际情况设置主频赫兹数
x_center = 0.5; % x 中心点
y_center = 0.5; % y 中心点
A = 1; % 振幅
att = 0.1; % 衰减系数

% temporal_temp_rk1
temporal_temp_rk1 = (1 - 2 * pi^2 * f_0^2 * t^2) * exp(-pi^2 * f_0^2 * t^2);

% temporal_temp_rk2
temporal_temp_rk2 = (1 - 2 * (pi * f_0 * (t - 0.2))^2) * exp(-(pi * f_0 * (t - 0.2))^2);

% temporal_temo_rk1d1t
temporal_temp_rk1d1t = (-4 * pi^2 * f_0^2 * t) * exp(-pi^2 * f_0^2 * t^2) + ...
                       (1 - 2 * pi^2 * f_0^2 * t^2) * exp(-pi^2 * f_0^2 * t^2) * (-2 * pi^2 * f_0^2 * t);

% temporal_temp_min
temporal_temp_min = t * exp(-3.5 * f_0 * t) * sin(2 * pi * f_0 * t);

% space_temp_Gauss
space_temp_Gauss = exp(-100 * ((x - x_center)^2 + (y - y_center)^2));

% space_temp_exp
space_temp_exp = A * exp(-att * sqrt((x - x_center)^2 + (y - y_center)^2));

r=temporal_temp_rk2 * space_temp_Gauss; 

%r=-sin(2*pi*x)*sin(2*pi*y)*exp(-t)